from pathlib import Path
import base64
import json
import math
import numpy as np
from PIL import Image, ImageDraw


OUT = Path(__file__).resolve().parent
PHOTO = Path(r"C:\Users\yopi1\Documents\01.Journals\08_EESD_Building_!!\Figures\EXP_Figures.jpeg")
COMPARISON = Path(r"C:\Users\yopi1\AppData\Local\Temp\codex-clipboard-8b182031-d296-4ee6-b404-3f7b6bcafd23.png")
RED = (255, 0, 0)

# Manually inspected source-pixel coordinates. Each visible section is a
# separate polyline; obscured portions are not extrapolated across steelwork.
TRACES = [
    ("left_column_lower_joint", [(418,1090),(429,1090),(439,1091),(451,1090)]),
    ("middle_infill_marked_joint", [(492,1025),(501,1024),(509,1023),(521,1023)]),
    ("lower_infill_upper_left", [(480,1092),(485,1092)]),
    ("lower_infill_upper_joint", [(521,1093),(534,1094),(543,1094),(552,1093),(563,1094),(580,1093)]),
    ("lower_infill_left_step", [(483,1187),(493,1188),(508,1189),(511,1192),(512,1197),(515,1202),(519,1204),(525,1203),(532,1202)]),
    ("lower_infill_main_step", [(545,1205),(552,1206),(556,1210),(561,1214),(562,1218),(570,1218),(575,1220),(580,1222),(583,1226),(583,1231),(586,1235),(591,1236),(598,1237),(600,1239),(599,1243)]),
    ("lower_infill_centre_rise", [(600,1238),(604,1231),(608,1227),(611,1225)]),
    ("lower_infill_right_step", [(621,1235),(620,1229),(620,1223),(627,1221),(635,1220),(640,1219),(646,1221),(658,1221),(672,1221)]),
    ("lower_infill_right_branch", [(639,1219),(640,1214),(641,1210)]),
    ("left_column_base_joint", [(419,1178),(429,1179),(433,1182),(444,1183),(450,1182)]),
    ("left_column_base_vertical", [(430,1185),(431,1194),(432,1202),(434,1209),(433,1215),(430,1221),(431,1232),(431,1237),(435,1240),(444,1241),(449,1240)]),
    ("left_column_open_joint_upper", [(418,1198),(425,1199),(429,1203)]),
    ("left_column_open_joint_lower", [(418,1219),(426,1218),(430,1218)]),
    ("left_column_foot_joint", [(418,1258),(429,1259),(439,1259),(446,1258)]),
    ("right_column_base_upper_joint", [(750,1216),(757,1216),(763,1217)]),
    ("right_column_base_upper_right", [(773,1217),(779,1218),(785,1219),(791,1221)]),
    ("right_column_open_step", [(765,1219),(763,1226),(761,1233),(761,1239),(750,1239),(750,1245),(753,1252),(756,1257),(754,1266),(755,1272)]),
    ("right_column_base_lower_joint", [(749,1272),(758,1273),(764,1274),(775,1275),(781,1276)]),
    ("right_column_base_cross_joint", [(761,1253),(769,1254),(779,1255),(788,1256)]),
]


def overlay(size, traces, width=2.2, transform=None):
    factor = 4
    layer = Image.new("RGBA", (size[0]*factor, size[1]*factor))
    draw = ImageDraw.Draw(layer)
    radius = width * factor / 2
    for _, points in traces:
        if transform:
            points = [transform(x, y) for x, y in points]
        points = [(x*factor, y*factor) for x, y in points]
        draw.line(points, fill=RED+(255,), width=round(width*factor), joint="curve")
        for x, y in points:
            draw.ellipse((x-radius,y-radius,x+radius,y+radius), fill=RED+(255,))
    return layer.resize(size, Image.Resampling.LANCZOS)


def fit_photo_registration(source, figure):
    # Fit the source photograph to its crop in the supplied comparison figure.
    # Coarse-to-fine color matching excludes text and the left numerical panel.
    a = np.asarray(source.convert("RGB"), dtype=float)
    b = np.asarray(figure.convert("RGB"), dtype=float)
    yy, xx = np.mgrid[55:558:9, 544:925:9]
    target = b[yy, xx]

    def error(sx_scale, sy_scale, tx, ty):
        sx = (xx-tx)/sx_scale
        sy = (yy-ty)/sy_scale
        xi = np.floor(sx).astype(int)
        yi = np.floor(sy).astype(int)
        if xi.min()<0 or yi.min()<0 or xi.max()>=a.shape[1]-1 or yi.max()>=a.shape[0]-1:
            return float("inf")
        dx = (sx-xi)[...,None]
        dy = (sy-yi)[...,None]
        sampled = (a[yi,xi]*(1-dx)*(1-dy) + a[yi,xi+1]*dx*(1-dy)
                   + a[yi+1,xi]*(1-dx)*dy + a[yi+1,xi+1]*dx*dy)
        return float(np.mean(np.minimum((sampled-target)**2, 3600)))

    best = (float("inf"), .496, 432., -154.)
    ranges = [(.006,4.,5.,.001,1.,1.), (.001,1.2,1.2,.0002,.2,.2), (.0002,.25,.25,.00005,.05,.05)]
    for rs,rx,ry,ds,dx,dy in ranges:
        _,s0,x0,y0=best
        candidates=[]
        for s in np.arange(s0-rs,s0+rs+ds/2,ds):
            for tx in np.arange(x0-rx,x0+rx+dx/2,dx):
                for ty in np.arange(y0-ry,y0+ry+dy/2,dy):
                    candidates.append((error(s,s,tx,ty),float(s),float(tx),float(ty)))
        best=min(candidates)
    # The supplied figure may have independently resized its photo width and
    # height. Refine all four affine parameters with a small Nelder-Mead fit.
    _,s,tx,ty=best
    initial=np.array([s,s,tx,ty])
    steps=np.array([.001,.001,1.,1.])
    simplex=np.vstack([initial]+[initial+np.eye(4)[i]*steps[i] for i in range(4)])
    values=np.array([error(*p) for p in simplex])
    for _ in range(240):
        order=np.argsort(values)
        simplex=simplex[order]
        values=values[order]
        if np.max(np.ptp(simplex,axis=0)/steps)<1e-5:
            break
        centroid=np.mean(simplex[:-1],axis=0)
        reflected=centroid+(centroid-simplex[-1])
        vr=error(*reflected)
        if values[0]<=vr<values[-2]:
            simplex[-1],values[-1]=reflected,vr
        elif vr<values[0]:
            expanded=centroid+2*(reflected-centroid)
            ve=error(*expanded)
            simplex[-1],values[-1]=(expanded,ve) if ve<vr else (reflected,vr)
        else:
            contracted=centroid+.5*((reflected if vr<values[-1] else simplex[-1])-centroid)
            vc=error(*contracted)
            if vc<min(vr,values[-1]):
                simplex[-1],values[-1]=contracted,vc
            else:
                simplex[1:]=simplex[0]+.5*(simplex[1:]-simplex[0])
                values[1:]=[error(*p) for p in simplex[1:]]
    i=np.argmin(values)
    return float(values[i]), *[float(v) for v in simplex[i]]


def save_svg(base_path, image_path, traces, stroke_width=2.2, transform=None):
    size = Image.open(image_path).size
    mime = "image/jpeg" if image_path.suffix.lower() in (".jpg", ".jpeg") else "image/png"
    encoded = base64.b64encode(image_path.read_bytes()).decode("ascii")
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{size[0]}" height="{size[1]}" viewBox="0 0 {size[0]} {size[1]}">',
             '<title>Experimental photograph with manually traced visible cracks</title>',
             '<desc>Red traces follow visible or previously marked cracks. Sections hidden by apparatus are left unconnected.</desc>',
             f'<image width="{size[0]}" height="{size[1]}" href="data:{mime};base64,{encoded}"/>',
             f'<g id="crack-traces" fill="none" stroke="#ff0000" stroke-width="{stroke_width}" stroke-linecap="round" stroke-linejoin="round">']
    for name, points in traces:
        if transform:
            points=[transform(x,y) for x,y in points]
        path="M "+" L ".join(f"{x:.3f} {y:.3f}" for x,y in points)
        lines.append(f'<path id="{name}" d="{path}"/>')
    lines.extend(["</g>","</svg>"])
    base_path.write_text("\n".join(lines), encoding="utf-8")


def main():
    photo = Image.open(PHOTO).convert("RGB")
    comparison = Image.open(COMPARISON).convert("RGB")
    photo_overlay = overlay(photo.size,TRACES)
    annotated = Image.alpha_composite(photo.convert("RGBA"),photo_overlay).convert("RGB")
    annotated.save(OUT/"EXP_Figures_red_cracks.png",dpi=(300,300))
    photo_overlay.save(OUT/"EXP_Figures_crack_overlay.png",dpi=(300,300))
    save_svg(OUT/"EXP_Figures_red_cracks_editable.svg",PHOTO,TRACES)

    loss,sx,sy,tx,ty=fit_photo_registration(photo,comparison)
    print(f"Registration: x_scale={sx:.6f}, y_scale={sy:.6f}, x_offset={tx:.3f}, y_offset={ty:.3f}, robust mean squared color error={loss:.3f}")
    transform=lambda x,y:(sx*x+tx,sy*y+ty)
    comp_overlay=overlay(comparison.size,TRACES,width=1.6,transform=transform)
    Image.alpha_composite(comparison.convert("RGBA"),comp_overlay).convert("RGB").save(OUT/"comparison_red_cracks.png",dpi=(300,300))
    save_svg(OUT/"comparison_red_cracks_editable.svg",COMPARISON,TRACES,stroke_width=1.6,transform=transform)

    annotated.crop((404,1037,800,1304)).resize((1188,801),Image.Resampling.LANCZOS).save(OUT/"crack_detail_preview.png")
    (OUT/"trace_coordinates.json").write_text(json.dumps({"source_size":photo.size,"source":str(PHOTO),"stroke":"#ff0000","source_stroke_width":2.2,"comparison_transform":{"sx":sx,"sy":sy,"tx":tx,"ty":ty},"traces":[{"id":n,"points":p} for n,p in TRACES]},indent=2))
    # Confirm the rendering preserves every source pixel outside the overlay.
    mask=np.asarray(photo_overlay.getchannel("A"))>0
    a=np.asarray(photo)
    b=np.asarray(annotated)
    assert np.array_equal(a[~mask],b[~mask])
    assert annotated.size==photo.size
    comparison_out=np.asarray(Image.open(OUT/"comparison_red_cracks.png"))
    comparison_mask=np.asarray(comp_overlay.getchannel("A"))>0
    assert np.array_equal(np.asarray(comparison)[~comparison_mask],comparison_out[~comparison_mask])
    assert np.array_equal(np.asarray(comparison)[:,:520],comparison_out[:,:520])
    print(f"Verified unchanged source dimensions and {int((~mask).sum()):,} pixels outside red traces.")


if __name__ == "__main__":
    main()

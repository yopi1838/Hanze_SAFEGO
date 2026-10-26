import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const pptxPath = "C:/Users/yopi1/Documents/09.Research_Dataset/Hanze_SAFEGO/00.Meeting/Progress_meeting_Route2_four_methods_2026-10-05_with_baseshear_and_method2_focus.pptx";
const outDir = "C:/Users/yopi1/Documents/09.Research_Dataset/Hanze_SAFEGO/.presentation_build_method2_focus/final-render";
await fs.mkdir(outDir, { recursive: true });
const p = await PresentationFile.importPptx(await FileBlob.load(pptxPath));
for (let i = 0; i < p.slides.items.length; i++) {
  const slide = p.slides.items[i];
  const png = await p.export({ slide, format: "png", scale: 1.5 });
  await fs.writeFile(path.join(outDir, `slide-${String(i + 1).padStart(2, "0")}.png`), new Uint8Array(await png.arrayBuffer()));
}
const montage = await p.export({ format: "webp", montage: true, scale: 1 });
await fs.writeFile(path.join(outDir, "montage.webp"), new Uint8Array(await montage.arrayBuffer()));
const inspection = await p.inspect({ kind: "slide,textbox,shape,image,chart,notes", maxChars: 30000 });
await fs.writeFile(path.join(outDir, "inspect.ndjson"), inspection.ndjson);
console.log(JSON.stringify({ slides: p.slides.items.length, outDir }, null, 2));

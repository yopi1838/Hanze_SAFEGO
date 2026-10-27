"""One figure for MECHANICAL_MODEL_US1.md: the model illustration (fig_mech_model_evolution.png, panels a-d) on top and
the pushover comparison (fig_dem_vs_mech.py in --combined mode, panels e-f) below, at the same width.
Run from mechanical_model_US1/."""
import os, sys, subprocess
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__))
py = sys.executable
subprocess.run([py, os.path.join(HERE, "fig_mech_model_evolution.py")], check=True, cwd=HERE)
top = Image.open(os.path.join(HERE, "fig_mech_model_evolution.png")).convert("RGB")
dpi = top.width / 13.5                                     # fig_dem_vs_mech is 13.5 in wide
subprocess.run([py, os.path.join(HERE, "fig_dem_vs_mech.py"), "--combined"], check=True, cwd=HERE,
               env=dict(os.environ, PO_DPI=str(dpi)), stdout=subprocess.DEVNULL)
bot = Image.open(os.path.join(HERE, "_pushover_panel.png")).convert("RGB")
if bot.width != top.width:
    bot = bot.resize((top.width, round(bot.height * top.width / bot.width)), Image.LANCZOS)
gap = 40
out = Image.new("RGB", (top.width, top.height + gap + bot.height), "white")
out.paste(top, (0, 0)); out.paste(bot, (0, top.height + gap))
out.save(os.path.join(HERE, "fig_mech_model_and_pushover.png"))
os.remove(os.path.join(HERE, "_pushover_panel.png"))
print("-> fig_mech_model_and_pushover.png", out.size)

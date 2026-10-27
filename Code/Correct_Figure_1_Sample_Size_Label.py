from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
import numpy as np
import cv2
import os

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "Figures/Source/Figure_1_author_supplied_original.png"
CORRECTED = ROOT / "Figures/Source/Figure_1_corrected_native.png"
FINAL = ROOT / "Figures/PNG/Figure_1.png"
FONT_CANDIDATES = [
    os.environ.get("FIGURE1_FONT_PATH", ""),
    "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf",
    "C:/Windows/Fonts/arial.ttf",
    "/Library/Fonts/Arial.ttf",
]
font_path = next((p for p in FONT_CANDIDATES if p and Path(p).is_file()), None)
if font_path is None:
    raise FileNotFoundError("Set FIGURE1_FONT_PATH to LiberationSans-Regular.ttf for the recorded output hash.")

image = Image.open(SOURCE).convert("RGB")
a = np.array(image)
mask = np.zeros(a.shape[:2], np.uint8)
x0, y0, x1, y1 = 35, 221, 290, 243
region = a[y0:y1, x0:x1]
gray = cv2.cvtColor(region, cv2.COLOR_RGB2GRAY)
local = (gray < 215).astype(np.uint8) * 255
local = cv2.dilate(local, np.ones((3, 3), np.uint8), iterations=1)
mask[y0:y1, x0:x1] = local
clean = cv2.inpaint(a, mask, 3, cv2.INPAINT_TELEA)
corrected = Image.fromarray(clean)
draw = ImageDraw.Draw(corrected)
font = ImageFont.truetype(font_path, 16)
draw.text((39, 222), "6 traits; 54,082 cases,", font=font, fill=(30, 34, 33))
corrected.save(CORRECTED, "PNG", dpi=(144, 144), optimize=True)
final = corrected.resize((round(corrected.width * 2.1), round(corrected.height * 2.1)), Image.Resampling.LANCZOS)
final.save(FINAL, "PNG", dpi=(300, 300), optimize=True)
print("Corrected only the combined disease-GWAS sample-size attribution; numerical totals were not changed.")

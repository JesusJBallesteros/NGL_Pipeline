"""Build a clean background plate: remove the bird + colored DLC dots from a photo
and inpaint the hole with surrounding floor texture (cv2 Telea).

    make_plate(photo_path, out_path, bbox=(x0,y0,x1,y1))

`bbox` is the bird's bounding box in photo pixels. Dark body OR colored dots inside
the box are masked, dilated to catch edges/shadow, then inpainted.
"""
import cv2
import numpy as np


def make_plate(photo_path, out_path, bbox, dark=85, sat=55, dilate=9, iters=2, radius=6):
    img = cv2.imread(photo_path)
    if img is None:
        raise FileNotFoundError(photo_path)
    H, W = img.shape[:2]
    x0, y0, x1, y1 = [int(v) for v in bbox]
    roi = img[y0:y1, x0:x1].astype(int)
    gray = roi.mean(axis=2)
    sat_roi = roi.max(axis=2) - roi.min(axis=2)
    m = (gray < dark) | (sat_roi > sat)
    mask = np.zeros((H, W), np.uint8)
    mask[y0:y1, x0:x1] = (m * 255).astype(np.uint8)
    k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (dilate, dilate))
    mask = cv2.dilate(mask, k, iterations=iters)
    plate = cv2.inpaint(img, mask, inpaintRadius=radius, flags=cv2.INPAINT_TELEA)
    cv2.imwrite(out_path, plate)
    return out_path


if __name__ == '__main__':
    import sys
    # usage: python make_plate.py photo.png out.png x0 y0 x1 y1
    p, o = sys.argv[1], sys.argv[2]
    box = tuple(map(int, sys.argv[3:7])) if len(sys.argv) >= 7 else (612, 158, 742, 266)
    print('wrote', make_plate(p, o, box))

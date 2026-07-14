import numpy as np

def normalize_angle(angle: float) -> float:
    return ((angle + 180.0) % 360.0) - 180.0

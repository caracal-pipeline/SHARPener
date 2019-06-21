import os
import sys

# get sharpener install directory
SHARPENER_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHARPENER_DIR = SHARPENER_PATH+'/sharpener/'

sys.path.append(os.path.join(SHARPENER_DIR, 'sharp_modules'))

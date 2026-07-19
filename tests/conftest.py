import sys, os

ROOT = os.path.dirname(os.path.abspath(__file__))
# make the app's scientific core importable
sys.path.insert(0, os.path.abspath(os.path.join(ROOT, "..", "app", "backend")))

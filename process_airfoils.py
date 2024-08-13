# %%

from airfoil import Airfoil
from pathlib import Path 

# dir_afs = Path("data/airfoils/djimatrice300rtk/")

fnames = dir_afs.glob("*.dat")

for fname in fnames:
    af = Airfoil.load_txt(fname)
    af.refine(101)
    af.normalize()
    af.save_txt(fname.with_suffix(".dat"))
    af.save_csv(fname.with_suffix(".csv"))
    af.plot(show=True)


# %%
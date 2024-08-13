
import AirfoilFast as AF

dir_out = "/home/hacs/projects/neaptide/Propeller/data/airfoils/dji9443"

names_afs = ["01", "02", "03", "04", "05", "06"]

for name_af in names_afs
    af = AF.Airfoil(joinpath(dir_out, name_af*".csv"))
    AF.save(af, joinpath(dir_out, name_af*".dat"))
end
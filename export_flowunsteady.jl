"""Convert Propeller.Rotor to FLOWUnsteady.Rotor"""

using Propeller
import YAML

# ==============================================================================
# Settings
# ==============================================================================

name_rotor = "dji_9443"
fname_config_rotor = "config/rotor/$name_rotor.yaml"
dir_out = "data/flowunsteady"
rpm = 1000

# ==============================================================================
# Convert DroneNoise.Rotor to FLOWUnsteady.Rotor
# ==============================================================================

config_rotor = YAML.load_file(fname_config_rotor)

rotor = Rotor(config_rotor)

export_flowunsteady(config_rotor, name_rotor, dir_out, rpm)


# Copy rotor config
# fname_out = "$dir_out/rotors/$name_rotor.yaml"
# run(`cp $fname_config_rotor $fname_out`)
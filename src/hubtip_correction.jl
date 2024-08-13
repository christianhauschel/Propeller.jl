
"""
Custom hub/tip-loss correction.

# Arguments 
- `r_tip`: tip radius [m]
- `r_hub`: hub radius [m]
- `r`: radial position [m]
- `inflow_angle`: inflow angle [rad]
- `n_blades`: number of blades [-]
- `t1`: tip correction parameter [-]
- `t2`: tip correction parameter [-]
- `t3`: tip correction parameter [-]
- `t_angle_min_d`: clipping threshold for the minimum allowable inflow angle [deg]
- `h1`: hub correction parameter [-]
- `h2`: hub correction parameter [-]
- `h3`: hub correction parameter [-]
- `h_angle_min_d`: clipping threshold for the minimum allowable inflow angle [deg]
"""
function hubtip_correction_custom(
    r_tip, r_hub, r, inflow_angle, n_blades;
    t1=1,
    t2=1,
    t3=1,
    t_angle_min_d=1,
    h1=1,
    h2=1,
    h3=1,
    h_angle_min_d=1
)
    inflow_angle_d = rad2deg(inflow_angle)

    if inflow_angle_d < t_angle_min_d
        F_tip = 1.0
    else
        f_tip = n_blades / 2 .* ((r_tip / r)^t1 - 1)^t2 / abs(sin(inflow_angle))^t3
        F_tip = 2 / π * acos(exp(-f_tip))
    end

    if inflow_angle_d < h_angle_min_d
        F_hub = 1.0
    else
        f_hub = n_blades / 2 .* ((r / r_hub)^h1 - 1)^h2 / abs(sin(inflow_angle))^h3
        F_hub = 2 / π * acos(exp(-f_hub))
    end
    return F_tip, F_hub
end

"""
Original Prandtl hub/tip-loss correction.
"""
function hubtip_correction_prandtl(r_tip, r_hub, r, inflow_angle, n_blades)
    return hubtip_correction_custom(
        r_tip, r_hub, r, inflow_angle, n_blades;
        t1=1,
        t2=1,
        t3=1,
        t_angle_min_d=1,
        h1=1,
        h2=1,
        h3=1,
        h_angle_min_d=1,
    )
end

"""
No hub/tip-loss correction.
"""
function hubtip_correction_none(r_tip, r_hub, r, inflow_angle, n_blades)
    return hubtip_correction_custom(
        r_tip, r_hub, r, inflow_angle, n_blades;
        t1=1,
        t2=0,
        t3=Inf,
        t_angle_min_d=5 * eps(),
        h1=1,
        h2=0,
        h3=Inf,
        h_angle_min_d=5 * eps(),
    )
end

"""
Modified Prandl hub/tip-loss correction with a strong hub correction.

t1 = 0.6
t2 = 5
t3 = 0.5
t_angle_min_d = 10

h1 = 2
h2 = 1
h3 = 0.25
h_angle_min_d = 0.05
"""
function hubtip_correction_modprandtl(r_tip, r_hub, r, inflow_angle, n_blades)
    return hubtip_correction_custom(
        r_tip, r_hub, r, inflow_angle, n_blades;
        t1=0.6,
        t2=5,
        t3=0.5,
        t_angle_min_d=10,
        h1=2,
        h2=1,
        h3=0.25,
        h_angle_min_d=0.05,
    )
end
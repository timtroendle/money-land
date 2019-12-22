import calliope
from calliope.core.util.logging import set_log_verbosity
import pyomo.core as po

ROOFTOP_TECH_NAME = "roof_mounted_pv"
UTILITY_TECH_NAME = "open_field_pv"
WIND_TECH_NAME1 = "wind_onshore_monopoly"
WIND_TECH_NAME2 = "wind_onshore_competing"
OFFSHORE_TECH_NAME = "wind_offshore"


def run(path_to_model, override_dict, roof_share, util_share, wind_share, path_to_output):
    assert roof_share + util_share + wind_share == 100
    set_log_verbosity("info", include_solver_output=True, capture_warnings=True)
    model = calliope.Model(
        path_to_model,
        scenario="no-hydro-costs,stylised-storage",
        override_dict=override_dict
    )
    model.run(build_only=True)
    pyomo_model = model.backend._backend
    pyomo_model.roof_constraint = po.Constraint(rule=rooftop_constraint(roof_share / 100))
    pyomo_model.util_constraint = po.Constraint(rule=utility_constraint(util_share / 100))
    pyomo_model.wind_constraint = po.Constraint(rule=wind_constraint(wind_share / 100))
    pyomo_model.offshore_constraint = po.Constraint(rule=offshore_constraint)
    results = model.backend.rerun().results
    results.attrs["scenario"] = f"roof-{roof_share}-percent,util-{util_share}-percent,wind-{wind_share}-percent"
    results.to_netcdf(path_to_output)


# FIXME make constraints location specific?!

def rooftop_constraint(share):
    def rooftop_constraint(model):
        lhs = sum(
            model.energy_cap[loc_tech]
            for loc_tech in model.loc_techs
            if ROOFTOP_TECH_NAME in str(loc_tech)
        )
        rhs = share * sum(
            model.energy_cap[loc_tech]
            for loc_tech in model.loc_techs
            if is_pv_or_wind(loc_tech)
        )
        return lhs == rhs
    return rooftop_constraint


def utility_constraint(share):
    def utility_constraint(model):
        lhs = sum(
            model.energy_cap[loc_tech]
            for loc_tech in model.loc_techs
            if UTILITY_TECH_NAME in str(loc_tech)
        )
        rhs = share * sum(
            model.energy_cap[loc_tech]
            for loc_tech in model.loc_techs
            if is_pv_or_wind(loc_tech)
        )
        return lhs == rhs
    return utility_constraint


def wind_constraint(share):
    def wind_constraint(model):
        lhs = sum(
            model.energy_cap[loc_tech]
            for loc_tech in model.loc_techs
            if is_wind(loc_tech)
        )
        rhs = share * sum(
            model.energy_cap[loc_tech]
            for loc_tech in model.loc_techs
            if is_pv_or_wind(loc_tech)
        )
        return lhs == rhs
    return wind_constraint


def offshore_constraint(model):
    lhs = sum(
        model.energy_cap[loc_tech]
        for loc_tech in model.loc_techs
        if OFFSHORE_TECH_NAME in str(loc_tech)
    )
    return lhs == 0


def is_wind(loc_tech):
    loc_tech = str(loc_tech)
    return (
        (WIND_TECH_NAME1 in loc_tech)
        or (WIND_TECH_NAME2 in loc_tech)
    )


def is_pv_or_wind(loc_tech):
    loc_tech = str(loc_tech)
    return (
        (ROOFTOP_TECH_NAME in loc_tech)
        or (UTILITY_TECH_NAME in loc_tech)
        or (WIND_TECH_NAME1 in loc_tech)
        or (WIND_TECH_NAME2 in loc_tech)
    )


if __name__ == "__main__":
    run(
        path_to_model=snakemake.input.model,
        override_dict=snakemake.params.override_dict,
        roof_share=int(snakemake.wildcards.roof),
        util_share=int(snakemake.wildcards.util),
        wind_share=int(snakemake.wildcards.wind),
        path_to_output=snakemake.output[0]
    )

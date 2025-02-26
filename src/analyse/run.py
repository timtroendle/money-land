import calliope
from calliope.core.util.logging import set_log_verbosity
import calliope.backend.run
import calliope.backend.pyomo
import calliope.core.attrdict
import calliope.exceptions
from calliope.analysis import postprocess
import pyomo.core as po

ROOFTOP_TECH_NAME1 = "roof_mounted_pv_n"
ROOFTOP_TECH_NAME2 = "roof_mounted_pv_e_w"
ROOFTOP_TECH_NAME3 = "roof_mounted_pv_s_flat"
UTILITY_TECH_NAME = "open_field_pv"
WIND_TECH_NAME1 = "wind_onshore_monopoly"
WIND_TECH_NAME2 = "wind_onshore_competing"
OFFSHORE_TECH_NAME = "wind_offshore"


def run(path_to_model, override_dict, roof_share, util_share, wind_share, offshore_share,
        units_without_shore, overrides, path_to_output):
    assert roof_share + util_share + wind_share + offshore_share == 100
    set_log_verbosity("info", include_solver_output=True, capture_warnings=True)
    model = calliope.Model(
        path_to_model,
        scenario=",".join(overrides),
        override_dict=override_dict
    )
    model.run(build_only=True)
    pyomo_model = model.backend._backend
    pyomo_model.roof_constraint = po.Constraint(pyomo_model.locs, rule=rooftop_constraint(roof_share / 100))
    pyomo_model.util_constraint = po.Constraint(pyomo_model.locs, rule=utility_constraint(util_share / 100))
    pyomo_model.wind_constraint = po.Constraint(
        pyomo_model.locs,
        rule=wind_constraint(wind_share / 100, offshore_share / 100, units_without_shore)
    )
    pyomo_model.offshore_constraint = po.Constraint(
        pyomo_model.locs,
        rule=offshore_constraint(offshore_share / 100, units_without_shore)
    )
    model = run_updated_model(model)
    scenario = f"roof-{roof_share}-percent,util-{util_share}-percent,wind-{wind_share}-percent,offshore-{offshore_share}-percent"
    model._model_data.attrs["scenario"] = scenario
    model.to_netcdf(path_to_output)


def rooftop_constraint(share):
    def rooftop_constraint(model, loc):
        lhs = sum(
            model.energy_cap[loc_tech]
            for loc_tech in model.loc_techs
            if is_rooftop_pv(loc_tech) and (loc_tech.split("::")[0] == loc)

        )
        rhs = share * sum(
            model.energy_cap[loc_tech]
            for loc_tech in model.loc_techs
            if is_pv_or_wind(loc_tech) and (loc_tech.split("::")[0] == loc)
        )
        return lhs == rhs
    return rooftop_constraint


def utility_constraint(share):
    def utility_constraint(model, loc):
        lhs = sum(
            model.energy_cap[loc_tech]
            for loc_tech in model.loc_techs
            if loc_tech.split("::") == [loc, UTILITY_TECH_NAME]
        )
        rhs = share * sum(
            model.energy_cap[loc_tech]
            for loc_tech in model.loc_techs
            if is_pv_or_wind(loc_tech) and (loc_tech.split("::")[0] == loc)
        )
        return lhs == rhs
    return utility_constraint


def wind_constraint(wind_share, offshore_share, units_without_shore):
    def wind_constraint(model, loc):
        if offshore_share > 0 and loc in units_without_shore:
            share = wind_share + offshore_share
        else:
            share = wind_share
        lhs = sum(
            model.energy_cap[loc_tech]
            for loc_tech in model.loc_techs
            if is_wind(loc_tech) and (loc_tech.split("::")[0] == loc)
        )
        rhs = share * sum(
            model.energy_cap[loc_tech]
            for loc_tech in model.loc_techs
            if is_pv_or_wind(loc_tech) and (loc_tech.split("::")[0] == loc)
        )
        return lhs == rhs
    return wind_constraint


def offshore_constraint(offshore_share, units_without_shore):
    def offshore_constraint(model, loc):
        if offshore_share > 0 and loc in units_without_shore:
            share = 0
        else:
            share = offshore_share
        lhs = sum(
            model.energy_cap[loc_tech]
            for loc_tech in model.loc_techs
            if loc_tech.split("::") == [loc, OFFSHORE_TECH_NAME]
        )
        rhs = share * sum(
            model.energy_cap[loc_tech]
            for loc_tech in model.loc_techs
            if is_pv_or_wind(loc_tech) and (loc_tech.split("::")[0] == loc)
        )
        return lhs == rhs
    return offshore_constraint


def is_wind(loc_tech):
    loc_tech = str(loc_tech)
    return (
        (WIND_TECH_NAME1 in loc_tech)
        or (WIND_TECH_NAME2 in loc_tech)
    )


def is_rooftop_pv(loc_tech):
    loc_tech = str(loc_tech)
    return (
        (ROOFTOP_TECH_NAME1 in loc_tech)
        or (ROOFTOP_TECH_NAME2 in loc_tech)
        or (ROOFTOP_TECH_NAME3 in loc_tech)
    )


def is_pv_or_wind(loc_tech):
    loc_tech = str(loc_tech)
    return (
        (ROOFTOP_TECH_NAME1 in loc_tech)
        or (ROOFTOP_TECH_NAME2 in loc_tech)
        or (ROOFTOP_TECH_NAME3 in loc_tech)
        or (UTILITY_TECH_NAME in loc_tech)
        or (WIND_TECH_NAME1 in loc_tech)
        or (WIND_TECH_NAME2 in loc_tech)
        or (OFFSHORE_TECH_NAME in loc_tech)
    )


def run_updated_model(model):
    # This method is largely taken from various places within Calliope's core code,
    # as Calliope does not offer this functionality.
    # The code is thus copyright Calliope authors.
    backend_model = model.backend._backend
    backend_model.__calliope_run_config = calliope.core.attrdict.AttrDict.from_yaml_string(
        model._model_data.attrs['run_config']
    )
    results, backend_mode = calliope.backend.run.run_plan(
        model_data=model._model_data,
        timings=model._timings,
        backend=calliope.backend.pyomo.model,
        backend_rerun=backend_model,
        build_only=False
    )

    # Add additional post-processed result variables to results
    if results.attrs.get('termination_condition', None) in ['optimal', 'feasible']:
        results = postprocess.postprocess_model_results(
            results, model._model_data, model._timings
        )
    else:
        raise calliope.exceptions.BackendError("Problem is non optimal.")

    for var in results.data_vars:
        results[var].attrs['is_result'] = 1
    model._model_data.update(results)
    model._model_data.attrs.update(results.attrs)
    model.results = model._model_data.filter_by_attrs(is_result=1)
    return model


if __name__ == "__main__":
    run(
        path_to_model=snakemake.input.model,
        override_dict=snakemake.params.override_dict,
        roof_share=int(snakemake.wildcards.roof),
        util_share=int(snakemake.wildcards.util),
        wind_share=int(snakemake.wildcards.wind),
        offshore_share=int(snakemake.wildcards.offshore),
        units_without_shore=snakemake.params.no_shore,
        overrides=snakemake.params.overrides,
        path_to_output=snakemake.output[0]
    )

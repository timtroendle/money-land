import pytest


EPSILON = 0.02 # 2% percent deviation is fine
ROOFTOP_TECHS = ["roof_mounted_pv_n", "roof_mounted_pv_e_w", "roof_mounted_pv_s_flat"]
UTILITY_TECH_NAME = "open_field_pv"
WIND_TECH_NAME1 = "wind_onshore_monopoly"
WIND_TECH_NAME2 = "wind_onshore_competing"
OFFSHORE_TECH_NAME = "wind_offshore"
WIND_TECHS = [WIND_TECH_NAME1, WIND_TECH_NAME2]
ALL_SUPPLY_TECHS = WIND_TECHS + ROOFTOP_TECHS + [OFFSHORE_TECH_NAME, UTILITY_TECH_NAME]
COUNTRIES_WITHOUT_SHORE = [ # FIXME should come from config
    "AUT",
    "BIH",
    "CHE",
    "CZE",
    "HUN",
    "LUX",
    "MKD",
    "SRB",
    "SVK",
    "SVN",
]


@pytest.fixture(scope="module")
def roof_share(model):
    return int(model.results.scenario.split(",")[0].split("-")[1]) / 100


@pytest.fixture(scope="module")
def util_share(model):
    return int(model.results.scenario.split(",")[1].split("-")[1]) / 100


@pytest.fixture(scope="module")
def wind_share(model):
    return int(model.results.scenario.split(",")[2].split("-")[1]) / 100


@pytest.fixture(scope="module")
def offshore_share(model):
    return int(model.results.scenario.split(",")[3].split("-")[1]) / 100


def test_roof_share(roof_share, energy_cap, location):
    roof_cap = energy_cap.sel(techs=ROOFTOP_TECHS, locs=location).sum("techs").item()
    all_cap = energy_cap.sel(techs=ALL_SUPPLY_TECHS, locs=location).sum("techs").item()
    assert roof_cap / all_cap == pytest.approx(roof_share, abs=EPSILON)


def test_util_share(util_share, energy_cap, location):
    util_cap = energy_cap.sel(techs=UTILITY_TECH_NAME, locs=location).item()
    all_cap = energy_cap.sel(techs=ALL_SUPPLY_TECHS, locs=location).sum("techs").item()
    assert util_cap / all_cap == pytest.approx(util_share, abs=EPSILON)


def test_wind_share(wind_share, offshore_share, energy_cap, location):
    if location in COUNTRIES_WITHOUT_SHORE:
        share = wind_share + offshore_share
    else:
        share = wind_share
    wind_cap = energy_cap.sel(techs=WIND_TECHS, locs=location).sum("techs").item()
    all_cap = energy_cap.sel(techs=ALL_SUPPLY_TECHS, locs=location).sum("techs").item()
    assert wind_cap / all_cap == pytest.approx(share, abs=EPSILON)


def test_offshore_share(offshore_share, energy_cap, location):
    if location in COUNTRIES_WITHOUT_SHORE:
        share = 0
    else:
        share = offshore_share
    offshore_cap = energy_cap.sel(techs=OFFSHORE_TECH_NAME, locs=location).item()
    all_cap = energy_cap.sel(techs=ALL_SUPPLY_TECHS, locs=location).sum("techs").item()
    assert offshore_cap / all_cap == pytest.approx(share, abs=EPSILON)

"""Create override requesting capacity shares."""
import jinja2

TEMPLATE = """
overrides:
    {%- for share in shares %}
    {%- for tech_name, tech_ids in techs.items() %}
    {{ tech_name }}-{{ (share * 100) | int }}-percent:
        group_constraints.{{ tech_name }}-{{ (share * 100) | int }}-percent:
            techs: {{ tech_ids }}
            demand_share_max:
                electricity: {{ share }}
    {%- endfor %}
    {%- endfor %}
"""
TECHS = {
    "wind": ["wind_onshore_monopoly", "wind_onshore_competing"],
    "util": ["open_field_pv"],
    "roof": ["roof_mounted_pv"]
}
SHARES = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] # FIXME inject


def capacity_share(path_to_result):
    """Generate a file that allows load shedding."""
    template = jinja2.Template(TEMPLATE)
    rendered = template.render(
        shares=SHARES,
        techs=TECHS
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    capacity_share(
        path_to_result=snakemake.output[0]
    )

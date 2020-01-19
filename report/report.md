# Introduction

Europe can generate its own electricity from renewable sources in Europe[@McKenna:2015; @Jacobson:2017; @Ruiz:2019; @Trondle:2019; @Bodis:2019a]. This provides an avenue to completely decarbonise the European electricity system. However, land use of renewable electricity is higher than land use of fossil and nuclear electricity[@Fthenakis:2009] (MORE SOURCES). The high land use may be problematic, first because it may lead to land use conflicts with other uses. Second, because generation infrastructure is perceived negatively by a part of the European population and high land uses of renewable electricity may be blocked by political barriers.

On the other hand, it has been shown that fully renewable electricity supply in Europe is possible in many different ways and that cost between very different designs must not vary much[@Schlachtberger:2017; @Neumann:2019; @Trondle:inreview]. While several studies have assessed cost[@Bussar:2014; @Connolly:2016; @Plessmann:2017; @Jacobson:2017; @Brown:2018; @Child:2019] and land use [@Fthenakis:2009; @Trondle:2019] of technologies and entire electricity systems, only one study has assessed the relationship between the two on the system level[@Palmer-Wilson:2019]. In their case study of Alberta, Canada, the authors find higher total system cost for lower land uses even though they allow for significant amounts of electricity from non-renewable sources. No study has assessed the relationship between land use and cost on the system level using only renewable resources in Europe. This leads me to the following research question: how can the land use of a fully renewable electricity system in Europe be minimised in the most cost-effective way?

Land use can be differentiated into exclusive land use and non-exclusive land use. For some technologies like onshore wind turbines, this distinction is important. The actual footprint of wind turbines is small compared to the spacing between them. The surface area between wind turbines can be used for other purposes, and is in fact often used for agriculture. Exclusive and non-exclusive land use of wind turbines are thus very different. Because I focus on political barriers based on public acceptance, instead of competition of land uses, I will consider the total land use in the following. For a discussion of the exclusive land use, see the Supplementary Material.

In the following, I will assess the land use and cost of fully renewable electricity systems in Europe with different shares of supply-side technologies: onshore wind, offshore wind, utility-scale PV, and rooftop PV. I use a nationally-resolved model of the electricity system to determine the necessary balancing capacities consisting of short-term (battery) and long-term (hydrogen) storage, bioenergy, and transmission. This allows me to determine the cost and land use for different shares of the four supply technologies.

# Results

## Land use of cost-optimal fully renewable electricity

Cost minimal system layout has 30% utility scale PV and 70% onshore wind -- no rooftop PV and no offshore wind. Such a system requires ~120,000 km^2^ or 2.5% of European land (Figure @fig:map). This is equivalent to the size of Switzerland and Austria.

![**Land use of solar and wind power in a cost-minimal, fully renewable electricity system with 30% utility-scale photovoltaics and 70% onshore wind capacity.** Countries considered in this study are shown in grey. The blue area shows the amount of land equivalent to the one that is required for solar and wind installations. In addition to solar and wind, electricity is generated from hydropower capacity of today's magnitude and bioenergy capacity generating electricity from biomass residuals.](report/map.png){#fig:map .class}

## Rooftop PV is not Pareto optimal in solar-wind systems

A system with only rooftop PV capacities requires no land at all, but its cost are 180% of the minimal cost (Figure {@fig:ternary-land-use-roof}a). More balanced combinations of onshore wind, utility-scale PV, and rooftop PV lead to cost and land-use combinations within this range: from the minimal cost up to 180% of the minimal cost, and from no land use to an area equivalent to the area of Switzerland and Austria. Only systems with no or almost no solar electricity have a land use even higher than that (Figure {@fig:ternary-land-use}b).

![**Cost and land use of fully renewable electricity systems with different shares of rooftop PV, utility-scale PV, onshore wind, and offshore wind.** **a,b,** Total system cost and land use of cases with only rooftop PV, utility-scale PV and onshore wind as supply-side options. **c,d** Total system cost and land use of cases with only offshore wind, utility-scale PV and onshore wind as supply-side options.](report/land-use/ternary.svg){#fig:ternary-land-use .class}

System designs including rooftop PV are almost never Pareto optimal; they are Pareto optimal only if the system is not at all supplied by wind (Figure @fig:scatter-land-use-roof). The land use of the cost optimal case can be reduced Pareto-optimally by exchanging wind capacity with utility-scale PV capacity. The first 10% land-use reductions are almost cost-neutral, and 50% land-use reductions can be achieved with a cost penalty of about 15%. Generally, the relationship is progressive, meaning that every additional land-use reduction is more costly than the one before. On average, reducing land use by moving from wind to utility-scale PV increases cost by 0.7 EUR / m^2^ / yr. Reducing land use further, by moving from utility-scale PV to rooftop PV, increases cost by 4.6 EUR / m^2^ / yr.

![**Cost and land use of fully renewable electricity systems with different shares of rooftop PV, utility-scale PV, and onshore wind.** All systems contain no offshore wind capacity. Cost is relative to minimal cost. Land use is relative to land use of the minimal-cost case. **a,** Cost and land use of all cases. The red line indicates the Pareto frontier on which neither cost nor land use can be reduced without increasing the respective other. **b,** Cost and land use of all cases. Darker colours indicate the share of rooftop PV (top), utility-scale PV (middle), and onshore wind (bottom).](report/land-use/scatter-roof.svg){#fig:scatter-land-use-roof .class}

## Offshore wind is cost-effective

Because rooftop PV is not Pareto optimal in solar-wind systems, I remove it and consider offshore wind as the third technology option in the following. Offshore wind is not included in a cost-minimal system design, but by enforcing high shares, onshore land use can be reduced to 10% of the cost-minimal case (Figure {@fig:ternary-land-use}d). Land use cannot be completely removed, because some countries have no shores and thus, they satisfy their domestic electricity demand with onshore wind instead of offshore wind (see Methods). Consequently, the European offshore share is never above 84%. In system designs with utility-scale PV and onshore and offshore wind, total system cost is never above 130% of the cost-minimal case, and cost is driven mainly by high shares of utility-scale PV (Figure {@fig:ternary-land-use}c).

Pareto optimal system designs always contain moderate amounts of utility-scale PV, and varying shares of onshore and offshore wind, leading to systems with either low cost or low land use (Figure @fig:scatter-land-use-offshore). Moving along the Pareto frontier from the cost-minimal case with 70% onshore wind to the land-use-minimal case with 70% offshore wind increases cost by 0.2 EUR / m^2^ / yr. Here, the curve is not progressive, but almost linear. The land-use-minimal case requires only about 12% of the onshore land of the cost-minimal case.

![**Cost and land use of fully renewable electricity systems with different shares of offshore wind, utility-scale PV, and onshore wind.** All systems contain no rooftop PV capacity. Cost is relative to minimal cost. Land use is relative to land use of the minimal-cost case. **a,** Cost and land use of all cases. The red line indicates the Pareto frontier on which neither cost nor land use can be reduced without increasing the respective other. **b,** Cost and land use of all cases. Darker colours indicate the share of offshore wind (top), utility-scale PV (middle), and onshore wind (bottom).](report/land-use/scatter-offshore.svg){#fig:scatter-land-use-offshore .class}

Onshore land use of fully renewable electricity systems can be most cost-effectively reduced using offshore wind. As long as offshore areas are not factored in, the cost is about 0.2 EUR / m^2^ / yr. The second most cost-effective option is utility-scale photovoltaics, which can reduce land use for about 0.7 EUR / m^2^ / yr. Rooftop PV, due to its high installation costs, is the least cost-effective option and requires up to 4.6 EUR / m^2^ / yr to reduce land use.

## Flexibility needs

Aside from different supply capacity, systems with different onshore land use require significantly different balancing capacity in terms of electricity storage, bioenergy, and transmission. Balancing needs are moderate for cases with balanced mixes of supply technologies (Figure @fig:flexibility). When supply is strongly biased towards one technology, flexibility needs rise, and in some cases they rise sharply. Exclusively-, or almost exclusively-solar cases require high amounts of short-term (battery) electricity storage. In the extreme case, storage capacities alone are able to fulfil the largest part of Europe's peak demand. Short-term storage capacities in these cases are combined with very high magnitudes of bioenergy capacity of up to 50% of peak demand to balance solar's seasonal fluctuations. Cases with mainly wind require much less bioenergy capacity and short-term storage capacity, but significantly more long-term storage capacities to balance wind fluctuations between days and weeks. In addition, they require around 2.5 times larger international transmission capacity than solar systems. While some of these numbers are very high, especially for cases with single supply technologies, there is no reason to believe these balancing capacities could not be built.

![**Flexibility needs for all analysed electricity system designs.** All designs are exclusively supplied from today's hydropower capacity and different shares of three additional technologies each: onshore wind and utility-scale PV in all cases, combined with either rooftop PV (top row) or offshore wind (bottom row). Each technology is varied from 0--100% of the total capacity of the three technologies. **a,b,c,d,** Storage power capacity (**a**), storage energy capacity (**b**), international transmission capacity (**c**), and bioenergy capacity (**d**). Not shown are hydropower capacities which are kept constant in all cases (36 GW run of river, 103 GW / 97 TWh reservoirs, 54 GW / 1.3 TWh pumped hydro storage).](report/flexibility.svg){#fig:flexibility .class}

Balancing technologies require land as well and thus add to the land use of the electricity system. Their land use is comparably small, however. [TODO Discuss land use of balancing using some number.]

## Uncertainty

Are the results robust towards uncertainties of technology cost and land use factors? I will perform a simple uncertainty analysis here. I will keep all system designs as they have been found, but I will vary the cost of supply technologies (and maybe balancing as well, but rather not) and their land use factors. This will allow me to show the impact of this uncertainty on the main results: 0.2 EUR / m^2^ / yr for offshore, 0.7 EUR / m^2^ / yr for utility-scale, and 4.6 EUR / m^2^ / yr for rooftop PV.

The approach is a simple approach, because the system design depends on the cost uncertainty. On the supply side, however, this relationship is not strong, because I enforce supply capacities. On the balancing side, this is more critical: reducing cost of bioenergy should lead to more bioenergy capacity, but using my approach, it won't. That's why I may not be able to analyse cost uncertainty of balancing technologies here.

I will, however, discuss the effect qualitatively. Also, I will discuss qualitatively the impact of electrifying the heat sector.

# Methods

I use the cost-minimising linear programming model introduced in ref [@Trondle:inreview]. I use a national spatial, and a 4h temporal resolution. I enforce national self-sufficiency, meaning that each country has to generate enough electricity to fulfil annual, domestic electricity demand.

Within this base model, for each considered case I enforce shares for four different supply technologies: onshore wind, offshore wind, utility-scale PV, and rooftop PV. I enforce relative installed capacity only, not the total amount of installed capacity. Because balancing losses are different in different cases, and because hydropower and biomass can supply electricity as well, the total installed capacity of wind and solar varies between the cases. For each case, the cost optimisation determines the cost-minimal total installed capacity of wind and solar, and the cost-minimal balancing capacities in each country. Balancing capacities comprise short-term storage (battery), long-term storage (hydrogen), bioenergy, and transmission.

I consider the land use factors as defined in ref [@Trondle:inreview] to determine land use based on installed capacity:

* onshore wind: 1 / 8 km^2^ / MW
* utility-scale PV: 1 / 80 km^2^ / MW
* offshore wind: 0 km^2^ / MW
* rooftop PV: 0 km^2^ / MW

# Bibliography

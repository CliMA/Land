# include("./scripts/check_weather_variation.jl")
weat_years = CSV.read("./data/Flagstaff/histo_30_amb.csv")
for year in 1971:2005
    weat_mask = (weat_years.Year .== year)
    weat_year = weat_years[weat_mask,:]
    mgst = mean(weat_year.Tair)
    sgsr = sum(weat_year.Rain)
    println(year, "\t", mgst, "\t", sgsr)
end

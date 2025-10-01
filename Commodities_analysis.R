library(readxl)

rm(list = ls())
data = read_excel("Dati_commodities.xlsx", col_types = c("text", rep("numeric", 18)))
data = as.data.frame(data)
rownames(data) = data[,1]
data = data[, -1]

# define meta-covariate
commodity_groups = list(
  Metals = c("Copper", "Platinum", "Silver", "Gold", "Aluminum", "Nickel"),
  Energy = c("OilPrice", "Coal"),
  Agri_Food = c("Soybean", "Cocoa", "Corn", "Onions", "Oranges", "Potatoes", "Rice"),
  Agri_Livestock_Fiber = c("Lambs", "Wool", "Cotton"))


map_to_group = function(x, groups) {
  for (g in names(groups)) {
    if (x %in% groups[[g]]) return(g)
  }
  return(NA)
}
commodity_category = sapply(colnames(data), map_to_group, groups = commodity_groups)
meta_covariate = data.frame(Commodity = colnames(data), Category = as.factor(commodity_category))
X = model.matrix(~meta_covariate$Category)

# maybe consider a subset of ten years or so after 1948 than there are not NA


### The analysis was performed with custom R (v. 3.5.3) code

### Packages used
library(ape);library(caper);
library(geiger);

### Pearson¡¯s correlations test
lm.function <- function(y,x,dt) {
  mm <- lm(dt[,y]~dt[,x],dt);
  a = as.numeric(format(coef(mm)[2]));
  b = as.numeric(format(coef(mm)[1]));
  r2 = as.numeric(format(summary(mm)$r.squared));
  pvalue<- as.numeric(cor.test(dt[,x], dt[,y], method = "pearson")$p.value);
  result.lm <- c(paste0("y= ", round(a,3),"x+",round(b,3)),
                 a, b, r2, pvalue)};

### Phylogenetic generalized least square model (PGLS); modified from paper doi 10.1007/s13127-015-0256-0
# Example:
# Relationship between the change in free energy of formation and the rate of thermal denaturation
dt_raw<- read.csv("./Input/data.csv", sep=",", header=T);
dt_raw$Energy<- -log10(-dt_raw$Energy);
dt_raw$lnReAct<- -(log10(-(dt_raw$lnReAct)));
dt_raw$Energy.lnReAct<- dt_raw$Energy/dt_raw$lnReAct;

Thermal.Index<- "lnReAct";
Trait<- "Energy";
data.new<- data.frame();
res.new<- data.frame();
tem<- 37;
dt<- subset(dt_raw, Temperature==tem);
tree<- read.tree("./Input/Tree.nwk");
row.names(dt)<- dt[,"species"];
check<- name.check(tree, dt); check;
tree<- drop.tip(tree, check$tree_not_data);
dt<- dt[!rownames(dt) %in% c(check$data_not_tree),];
sum(sort(tree$tip.label) != sort(dt$species));
tree$node.label<- NULL;   # Removing node labels
  
comp.data<- comparative.data(phy=tree, data=dt, names.col="species");
mod.pgls<- pgls(lnReAct ~ Energy, data=comp.data);
aic<- mod.pgls$aic;
a<- as.numeric(format(coef(mod.pgls)[2], digits = 8));
b<- as.numeric(format(coef(mod.pgls)[1], digits = 8));
res.func<- paste0("y= ", round(a,3),"x+",round(b,3));
pvalue<- as.numeric(format(1- pf(summary(mod.pgls)$fstatistic[1], summary(mod.pgls)$fstatistic[2], summary(mod.pgls)$fstatistic[3]), digits = 8));
r2<- as.numeric(format(summary(mod.pgls)$r.squared, digits = 3));
res.func<- paste0("y= ", round(a,3),"x+",round(b,3));

res.new<- cbind.data.frame(res.func,a,b,round(r2,3),round(pvalue,3),tem,Thermal.Index,Trait,aic,paste0(Thermal.Index,"~",Trait));   # Model results
data.new<- cbind.data.frame(mod.pgls$fitted,(mod.pgls$fitted-b)/a,
                            mod.pgls$data$data$lnReAct,mod.pgls$data$data$Energy,
                            tem,
                            paste0(Thermal.Index,"~",Trait));   # Fitted values
data.new$species<- row.names(data.new);

names(data.new)<- c(paste0(Thermal.Index,"_pgls"),paste0(Trait,"_pgls"),Thermal.Index,Trait,"Temperature","model.Call","species");
colnames(res.new)<- c("function","a","b","R2","pValue","Temperature","Thermal.Index","Trait","AIC","model.Call");
res.new<- res.new[order(res.new[,"Temperature"],decreasing=F),];



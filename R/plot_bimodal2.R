library(ggplot2)
library(latex2exp)

wd_vec <- c("simResults",
            "simResultsrhopos05",
            "simResultsrhopos08",
            "simResultsrhoneg05",
            "simResultsrhoneg08")
for(wd in wd_vec){
setwd(paste("<LOCAL LOCATION OF ZDIRECT>",wd, sep = ""))

method_choice = c(
  "zdirect",
  "Storey",
  "StoreyAdaptive",
  "BH",
  "ash",
  "lfsr",
  "GR",
  "dBH")
method_names = c(
  "zdirect",
  "Storey",
  "StoreyAdaptive",
  "BH",
  "ash",
  "lfsr",
  "GR",
  "dBH")
method_labels = c(
  "ZDIRECT",
  "Storey_dir",
  "aStorey_dir",
  "BH_dir",
  "ASH",
  "LFSR",
  "GR",
  "dBH")
method_colors = c(
  "red",
  "purple", "#00BCE3",
  "blue", "lightgreen",  "plum4", "orange", "hotpink"
)

names(method_colors) = method_names

method_labels = unname(TeX(c("ZDIRECT", "$STS_{dir}$", "$aSTS_{dir}$", "$BH_{dir}$", "ASH", "LFSR", "GR", "dBH" )))

scale_shape_values = 0:(length(method_names)-1)
scale_linetype_values = c(1,1,1, 2:(length(method_names)-2))

nullprop_choice = c(0, 0.2,  0.5, 0.8)
symm_choice = c(0.5, 0.75, 1)
# nullprop_choice = c( 0.8)
# symm_choice = c(0.9)

mu_choice = seq(0.5, 2.5, by = 0.5)

nullprop_vec = c()
mu_vec = c()
symm_vec = c()
method_vec = c()
FDR = c()
TPR = c()


for (nullprop in nullprop_choice){
  for(mu in mu_choice){
    for (symm in symm_choice){
      for (method in method_choice){
        filepath = paste(
          "/nullprop", nullprop,
          "mu", mu,
          "symm", symm ,
          ".RData", sep = "" )
        load(paste0( getwd(),  filepath))
        nullprop_vec  = c(nullprop_vec , nullprop)
        mu_vec  = c(mu_vec , mu)
        symm_vec = c(symm_vec, symm)
        method_vec = c(method_vec, method)
        FDR = c( FDR,  mean(eval(as.symbol(paste0("FDP_",  method)))))
        TPR = c( TPR,  mean(eval(as.symbol(paste0("TPP_",  method)))))
      }
    }
  }
}

data = data.frame(nullprop_vec, mu_vec, symm_vec, method_vec, FDR, TPR)

# symm_names= paste0("symm = ", symm_choice)
symm_vec <- as.factor(symm_vec)
data$symm_vec <- as.factor(data$symm_vec)
nullprop_vec <- as.factor(nullprop_vec)
data$nullprop_vec <- as.factor(data$nullprop_vec)

symm_names_list <- list(
  '0.5'=TeX(c("$v = 0.5$")),
  '0.75'=TeX(c("$v = 0.75$")),
  '1'=TeX(c("$v = 1$"))
)

nullprop_names_list <- list(
  '0'=TeX(c("$w = 0$")),
  '0.2'=TeX(c("$w = 0.2$")),
  '0.5'=TeX(c("$w = 0.5$")),
  '0.8'=TeX(c("$w = 0.8$"))
)

symmnullprop_labeller <- function(variable,value){
  if (variable=='symm_vec') {
    return(symm_names_list[value])
  } else {
    return(nullprop_names_list[value])
  }
}


grid_2d_FDR = facet_grid( nullprop_vec ~ symm_vec,
                          scales = "free",
                          labeller = symmnullprop_labeller)

grid_2d_TPR = facet_grid(nullprop_vec ~ symm_vec,
                         scales = "free",
                         labeller = symmnullprop_labeller)


#
# ggplot_list = list()
# count =0
# for (target in  c("FDR", "TPR")){
#   count = count +1
#   if (target == "FDR") {
#     # ylab_name = "False discovery rate"
#     ylab_name = ""
#     out_name = paste0("~/Dropbox/punim1304/ZDIRECT/icml2022/FDR_Plot.pdf")
#     legend_pos = "bottom"
#   }
#   if (target == "TPR") {
#     # ylab_name = "True positive rate"
#     ylab_name = ""
#     out_name = paste0("~/Dropbox/punim1304/ZDIRECT/icml2022/TPR_Plot.pdf")
#     legend_pos = "bottom"
#   }



ggobj_FDR = ggplot(data,
                   aes(x=mu_vec,
                       y = FDR,
                       group=method_vec,
                       colour = method_vec,
                       shape = method_vec)) +
  geom_line(aes(linetype=method_vec))+
  geom_point(aes(shape=method_vec)) +
  xlab(expression(xi)) +
  ylab("Directional false discovery rate") +
  scale_color_manual(name="",
                     breaks = method_names,
                     labels= method_labels,
                     values=method_colors) +
  scale_shape_manual(name="",
                     breaks = method_names,
                     labels= method_labels,
                     values=scale_shape_values)+
  scale_linetype_manual(name="",
                        breaks = method_names,
                        labels=method_labels,
                        values=scale_linetype_values)+
  geom_hline(yintercept=0.10) +
  ylim(c(0,0.20)) +
  grid_2d_FDR+
  guides(colour = guide_legend(nrow = 1)) +
  theme_bw()+
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))




ggobj_TPR = ggplot(data,
                   aes(x=mu_vec,
                       y = TPR,
                       group=method_vec,
                       colour = method_vec,
                       shape = method_vec)) +
  geom_line(aes(linetype=method_vec))+
  geom_point(aes(shape=method_vec)) +
  xlab(expression(xi)) +
  ylab("True positive rate") +
  scale_color_manual(name="",
                     breaks = method_names,
                     labels= method_labels,
                     values=method_colors) +
  scale_shape_manual(name="",
                     breaks = method_names,
                     labels= method_labels,
                     values=scale_shape_values)+
  scale_linetype_manual(name="",
                        breaks = method_names,
                        labels=method_labels,
                        values=scale_linetype_values)+
  grid_2d_TPR +
  guides(colour = guide_legend(nrow = 1)) +
  theme_bw()+
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))






ggobj = ggpubr::ggarrange(ggobj_FDR, ggobj_TPR, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")

ggobj

# if (target == "TPR"){
#   fig_height = 4
#   fig_width = 4
# }
# if (target == "FDR"){
#   fig_height = 4
#   fig_width = 4
# }

ggobj
# save picture
out_name = paste0(wd,".pdf",collapse = "")
pdf(file = out_name,
    width = 8.5, height = 11.5)
ggobj
dev.off()

}

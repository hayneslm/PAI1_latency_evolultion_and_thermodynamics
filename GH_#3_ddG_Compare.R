# This code compared ddG values calculated using differnt methods and PAI- structures

setwd("~/Documents/R_projects/Binning_2024-06/")

nbins = 50

# Import active state ddG data for EvoEF and xFold
ddG_3q02 = read_tsv("new_ddG/active_3Q02_dG_value.xls") %>% 
  select(Ref_AA, Pos, Mut_AA, Evo_REF, Evo_Ref_WT, Evo_Mut, 
         foldx_Ref_WT, foldx_Mut) %>% 
  rename(AA = Ref_AA, Mut = Mut_AA) %>% 
  mutate(Variant = paste0(AA, Pos, Mut), 
         EvoEF = Evo_Mut-Evo_Ref_WT,
         xfold = foldx_Mut-foldx_Ref_WT) %>% 
  select(AA, Pos, Mut, Variant, EvoEF, xfold)

model <- lm(xfold ~ EvoEF, data = ddG_3q02)
summary(model)

ggplot(ddG_3q02, aes(x = EvoEF, y = xfold))+
  geom_jitter(color = 'grey60')+
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed")+
  theme_classic()

# Write a pdf of active EvoEF vs xFold
ggsave("Active_EvoEF_vs_xFold.pdf", device = "pdf")

# Import latent state ddG data for EvoEF and xFold
ddG_1dvn = read_tsv("new_ddG/latent_1DVN_dG_value.xls") %>% 
  select(Ref_AA, Pos, Mut_AA, Evo_REF, Evo_Ref_WT, Evo_Mut, 
         foldx_Ref_WT, foldx_Mut) %>% 
  rename(AA = Ref_AA, Mut = Mut_AA) %>% 
  mutate(Variant = paste0(AA, Pos, Mut), 
         EvoEF = Evo_Mut-Evo_Ref_WT,
         xfold = foldx_Mut-foldx_Ref_WT) %>% 
  select(AA, Pos, Mut, Variant, EvoEF, xfold)



model <- lm(xfold ~ EvoEF, data = ddG_1dvn)
summary(model)

ggplot(ddG_1dvn, aes(x = EvoEF, y = xfold))+
  geom_jitter(color = 'grey60')+
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed")+
  theme_classic()

# Write a pdf of latent EvoEF vs xFold
ggsave("Latent_EvoEF_vs_xFold.pdf", device = "pdf")

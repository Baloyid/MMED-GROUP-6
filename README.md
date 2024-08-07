# How do agent based models and compartmental models differ in their predictions of SIRS transmission dynamics ?

Overall research topic: Investigating model assumptions using model comparisons. 

Members: Dumisane, Lindokuhle, Yanchao, Ryan, Wisdom and Stanley. 

Facilitator: Slimane BenMiled

---

![SIRS Model Diagram](SIRSModel.png)

---

The repository contains the code and work for group 6. The repository is structured into folders for code, data, outputs and resources.


1. The **data** folder contains folders holding data for different simulations titled `pX`, where X is the population size for the datasets contained within it, so that, p500 implies all datasets within the p500 folder are from ABM simulations where the population size is fixed at 500. Every folder within the data folder also has a README which contains details of the ABM simulations. Within each data folder are three more folders: a no-neighbours, four-neighbours and eight-neighbours folder. Reeaders should consult the img/neighbourhood.png image to have a look at what the neighbourhood structures are. Thus for each population size we have a total of 30 datasets.
    
2. The **code-files** folder contains R and netlogo code files. There are seperate files for ODE-SIRS and ChainBinomial-SIRS. The files are well structured to ensure that the code flow is intuitive, and systematic. In every code file (sirs-ode, sirs-chainbinomial), the resulting output is stored in the outputs folder
   
3. The **outputs** folder contains the outputs from running the SIRS models (ODE & ChainBinomial). The output has seperate folders for the different models used. The csv files labelled using similar method of `pX`.

4. The **img** folder contains images that were generated while running the analysis, as well as model diagrams.

5. The **globalenv** folders contain a snapshot rdata for a clean run of the ODE and chainbinomial runs. Loading this into a fresh R session gives you all the data that was produced when running the optimizations.

---

To be able to reproduce the analysis results and figures, one needs to either run the sirs-ode.R file to get ODE results or sirs-chainbinomial.R to get the chain-binomial results.




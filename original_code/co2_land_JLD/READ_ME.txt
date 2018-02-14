READ_ME: LandUptakeModel_Driver.m 
(formerly “nonlin_land_Qs_annotate.m”)
by Lauren Rafelski, lauren.rafelski@gmail.com, Oct 25, 2012, modified by Julia Dohner beginning September 2017

LandUptakeModel_Driver.m is the driver code for the land model. This code fits a model of land uptake to the residual land flux. The model includes CO2 fertilization or N fertilization, and temperature-dependent respiration.

From Rafelski (2009) pg. 721:
“Q10 is the factor that the fluxes change for a 10◦C temperature change.”

We also consider simpler cases in which one or more processes are neglected, for example, by setting epsilon = 0, gamma = 0, or Q10 = 1.

The two-box model that we use is a more mechanistic approach to understanding variability, allowing for time constants in the response of CO2 fluxes to forcings. The parameters epsilon, gamma and Q10 are fitted empirically using observed atmospheric CO2 and carbon emissions.”

Temperature-independent model: “We now compare estimates of the residual flux to results from a simple land biosphere model, which consists of the biospheric model driven only by CO2 fertilization (gamma = 1 and Q10 = 1).

Temperature-dependent model: “We now repeat the same analysis, allowing for temperature-dependent respiration (relaxing the requirement that Q10 = 1) driven by changes in the global land air temperature record. We thus now fit both epsilon and Q10.”

CO2 or N fertilization: Q10 either defined as 1 or fit to land air temp record

Some notes:
Lines 16-23: define what kind of model fit you want to do: High or low land use emissions, nitrogen fertilization or not, 10-year filter on the data or not

Lines 26-30: load different temperature records (you might want updated records)

Lines 40-50: CO2 fluxes, atmospheric CO2 (you might want to update)

Line 81: moving boxcar average of the land temperature. Note: this was from earlier tests, and this land record doesn't end up being used in the final runs

Lines 98-131: chooses a temperature record to use. Global land air temperature was the basic case; also tried NPP-weighted temperatures

Lines 158 and 162: nlinfit is a routine that comes with MATLAB. If you don't have it, let me (LR) know

		land_fit_Qs_annotate.m: you have to go into this subroutine to change 			between temperature-dependent and temperature-independent runs, or between 		CO2 fertilization and N fertilization runs.

		bioboxtwo_sub10_annotate.m and bioboxtwo_subN_annotate: Choose between 			lines 44 and 46 and lines 50 and 52 to switch between temperature-			dependent respiration and temperature-dependent photosynthesis

Lines 166-176: calculates covariance and correlation. You may need to change the index numbers for J and resid on line 169, depending on the range you want to look at

Lines 181-200: the model is run again with the best fit values so you can plot or save the results. Choose between the CO2 model on line 196 or the N2 model on line 200 (you could also make this automatic based on the parameter value of "nitrogen" (line 21)
# Leave-one-out-cross-validation-of-EBV
This Julia script was developed to implemment an efficient LOOCV method for prediction of breeding values and other random effects 
under a general mixed linear model with multiple random effects. Details can be seen:
https://onlinelibrary.wiley.com/doi/full/10.1111/jbg.12545 

To test the LOOCV function (LOOCV_example.jl), use the small example data from the data folder. 
The LOOCV test run is based on the model: Nur2ADG= Batch + Died + EntryAge + NurPenBatch + SowID + e, where Batch and Died are the fixed effect, EntryAge is the covariate, NurPenBatch and SowID are the random effect.

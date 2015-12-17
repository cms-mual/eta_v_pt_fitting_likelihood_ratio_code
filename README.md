# eta_v_pt_fitting_likelihood_ratio_code
Code to test the likelihood ratio of a three term fit vs a 2 term fit for our eta distributions of pt slices

just run python likelihood_ratio_implementation.py -b to replicate my results

there are lots of paramaters to play around with in the code (like the number of bins) 

Also, the likelihood fitting code is a bit finiky. Check the Chi2/ndf or plot to visualy make sure it is fitting as expected.

Read more about the math here: 

https://en.wikipedia.org/wiki/Likelihood-ratio_test

http://ac.els-cdn.com/0167508784900164/1-s2.0-0167508784900164-main.pdf?_tid=6b4698ba-a415-11e5-b84e-00000aab0f6c&acdnat=1450284940_492d28bf72d4b5f92541f61ac179fd37

https://en.wikipedia.org/wiki/Chi-squared_distribution

http://maxwell.ucsc.edu/~drip/133/ch4.pdf


warning, lots of plots will be produced! Feel free to comment out the draw or save as parts




NOTES: 

Structural mistakes in this lib:

0. The way the arrays are fed to and output by scipy.optimize.leastsq are ridiculous. I should have extracted the data by name e.g. in a pnadas.Series as soon as I received the output of leastsq. I did not, which resulted in code I know I won't understand a few months from now for extraction of data from those arrays ( specifically pf, pferr ).

1. It was unwise to embed the natural 2D structure of the data ( x and y pixel coords for each configuration ) in a 1D table, with separate cols for x and y. This made life harder later on. Also, since there are multiple files with the same DB structure, the dataset as a whole actually has a 3D structure and would naturally all fit in the same DB. I regret not doing this.

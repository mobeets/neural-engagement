Code for analyses and figures for paper ["Learning is shaped by abrupt changes in neural engagement"](https://www.nature.com/articles/s41593-021-00822-8) by Hennig et al. (2021)

## Requirements

This code was developed using Matlab R2015a, along with the [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html) and the [Optimization Toolbox](https://www.mathworks.com/products/optimization.html). The code has also been tested on Matlab R2019a, R2019b, and R2020a. Code may not be compatible with Matlab versions earlier than R2014b.

## Usage

This codepack includes data from two example sessions (20120528 from monkey J, and 20160728 from monkey N; see `data/preprocessed`). This code will generate data figures similar to those in the main text, but using the provided example sessions rather than all sessions as in the paper. To generate these figures, simply run the following:

```matlab
>> figure.F2; % Fig. 2
>> figure.F3; % Fig. 3
>> figure.F4; % Fig. 4
>> figure.F5and6; % Fig. 5 and Fig. 6
```

The above code takes around 20 seconds to run on a Macbook Pro.

## Instructions for use

© Jay Hennig, 2020

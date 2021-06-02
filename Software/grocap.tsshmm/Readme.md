FivePrime is a method devloped by the Siepel and Lis labs at Cornell University. Here we only included codes that we modified, in order to run our pipeline on multiple methods. The original code can be retrieved from [GitHub](https://github.com/andrelmartins/grocap.tsshmm).

Modified files:
* `rpkg/grocaptss/pkg/R/hmm.parse.R`: Added fail-safe strategy for cases where viterbi failed;
* `rpkg/grocaptss/pkg/R/postproc.R`: Fixed short-circuited `thresh` in `create.pairs()`.

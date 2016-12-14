# Small framework to derive L2 residual corrections

First derive the response using deriveResponse.C
Then sum the histos and do some plots using sumResponse.C

To check closure, use deriveClosure.C

All of the above is possibly to launch in batch mode on lxplus. I used a workaround in order to control 
the submission with python. It basically generates separate file lists and .csh scripts for submission.
Then it creates directories on eos to store the output.

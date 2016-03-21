# batch_experiment

Things you could want to do, and this tool does for you:
* You want to run a batch of sh commands, but only one of those per core/cpu.
* You want to maximize the core use, the moment a core/cpu becomes free from one of your commands, you want the next command to take its place.
* You want the output of those commands to be saved by command/file (you want to have a file with what exactly the command outputted).
* You want to specify timeouts for those commands to be killed.
* You want the power to resume the batch from an interruption (i.e. system crashes, or energy goes down) without losing any work.

What conditions you need to use this tool:
* You use linux.
* You have installed: sh (the shell); time (NOT the bash/sh internal command, but the package one, i.e.: https://www.archlinux.org/packages/extra/x86_64/time/ ); timeout (from coreutils); taskset (from util-linux, https://www.archlinux.org/packages/core/x86_64/util-linux/ ).

What is not needed:
* To know how to program in ruby. Only taking less than 5 minutes to learn some basic syntax will suffice. However, there's an exception, if you want not only to execute the commands but want to extract and group some information from the output to a CVS, you will need to tell ruby how to do the extracting part.

## Example



This code was born in [this repository](https://github.com/henriquebecker91/masters/tree/master/codes/rb/batch_experiment).


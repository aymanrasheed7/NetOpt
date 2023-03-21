# NetOpt
## General command
<pre>
outFile ARG1 ARG2 ARG3 OPTIONALARG4 OPTIONALARG5
</pre>
## Argument examples
<pre>
ARG1: "E12" (using {1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2}), "ONE", "INT", "ODD"
ARG2: "E", "PI", "PHI", "SQRT" (sqrt(nR)), "SQRT2.5", "2.5", "-2.5" (2.5 with optional exclution)
ARG3: "1", "2", ..., "12" (nR = number of resistors)
ARG4: "1", "2", ..., "12" (table size controller)
ARG5: "1", "2", ..., "256" (number of threads)
</pre>
## Example usage
### Windows
<pre>
g++ .\NetOptThread.cpp; .\a.exe E12 SQRT 10
</pre>
### Linux
<pre>
g++ -pthread NetOptThread.cpp; ./a.out E12 SQRT 10
</pre>
## Example output
### Windows
<pre>
Solution: (1.2+(3.9|(2.2+(((1.0|(1.5+4.7))|(2.7|5.6))+(1.8|3.3)))))
Result: 3.1622776652262714 (845658744/267420775)
Target: 3.1622776601683795
Cost: 5.057892e-09
CPU time (s): 12.819000000
Execution time (s): 12.815042700
</pre>
### Linux
<pre>
Solution: (1.2+(3.9|(2.2+(((1.0|(1.5+4.7))|(2.7|5.6))+(1.8|3.3)))))
Result: 3.1622776652262714 (845658744/267420775)
Target: 3.1622776601683795
Cost: 5.057892e-09
CPU time (s): 135.711018000
Execution time (s): 8.269306404
</pre>

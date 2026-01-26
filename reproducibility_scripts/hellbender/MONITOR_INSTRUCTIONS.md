# How to Monitor SKIOME Job in Real-Time

## Option 1: Use the Monitoring Script (Easiest)

Run this in your terminal:

```bash
cd /Users/nathanbresette/Desktop/MeLSI/reproducibility_scripts/hellbender
bash monitor_skiome.sh
```

This will show live output from the job. Press `Ctrl+C` to stop.

## Option 2: Manual SSH with tail -f

Run this in your terminal:

```bash
ssh nbhtd@hellbender.rnet.missouri.edu "cd ~/melsi_simulations/hellbender && tail -f skiome_run2.err"
```

Or to see both output and error:

```bash
ssh nbhtd@hellbender.rnet.missouri.edu "cd ~/melsi_simulations/hellbender && tail -f skiome_run2.out skiome_run2.err"
```

## Option 3: Check Progress Periodically

```bash
ssh nbhtd@hellbender.rnet.missouri.edu "cd ~/melsi_simulations/hellbender && tail -20 skiome_run2.err"
```

## What to Look For

- **"Generating null distribution with 200 permutations..."** - Running permutations (slow step)
- **"[Permutation X of 200]"** - Progress through permutations
- **"MeLSI Results:"** - Analysis complete
- **"Results saved to: skiome_validation_results.csv"** - Job finished!

## Check if Job is Still Running

```bash
ssh nbhtd@hellbender.rnet.missouri.edu "pgrep -f skiome && echo 'Job is running' || echo 'Job is not running'"
```

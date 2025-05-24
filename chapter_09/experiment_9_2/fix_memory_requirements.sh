#!/bin/bash

echo "🔧 Fixing Nextflow memory requirements for low-memory systems..."

# Update base.config with lower memory requirements
cat > conf/base.config << 'BASEEOF'
/*
 * -------------------------------------------------
 *  Base Nextflow config file for low-memory systems
 * -------------------------------------------------
 */

process {
    cpus   = 1
    memory = 2.GB
    time   = 2.h

    errorStrategy = { task.exitStatus in [143,137,104,134,139] ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'

    withLabel:process_single {
        cpus   = 1
        memory = 1.GB
        time   = 1.h
    }
    withLabel:process_low {
        cpus   = 1
        memory = 2.GB
        time   = 2.h
    }
    withLabel:process_medium {
        cpus   = 2
        memory = 3.GB
        time   = 4.h
    }
    withLabel:process_high {
        cpus   = 2
        memory = 4.GB
        time   = 6.h
    }

    // Specific process configurations with low memory
    withName:RUST_PSEUDOALIGN {
        cpus   = 2
        memory = 3.GB
        time   = 4.h
    }

    withName:CUSTOM_DUMPSOFTWAREVERSIONS {
        cpus   = 1
        memory = 1.GB
        time   = 1.h
    }
}
BASEEOF

echo "✅ Memory requirements fixed!"
echo "Now run: nextflow run main.nf --method pseudo --reads 'data/*.fastq*' -resume"

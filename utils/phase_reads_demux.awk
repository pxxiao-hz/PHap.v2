BEGIN {
    FS = OFS = "\t"
    for (quality_code = 33; quality_code <= 126; quality_code++) {
        phred_by_character[sprintf("%c", quality_code)] = quality_code - 33
    }

    while ((getline line < groups_file) > 0) {
        count = split(line, fields, "\t")
        if (count != 2 || fields[1] == "" || fields[2] == "") {
            print "PHAP_ERROR\tinvalid group FIFO row: " line > "/dev/stderr"
            exit 2
        }
        fifo_by_group[fields[1]] = fields[2]
    }
    close(groups_file)
    if (length(fifo_by_group) == 0) {
        print "PHAP_ERROR\tno output groups" > "/dev/stderr"
        exit 2
    }

    # Opening every FIFO here guarantees that every compressor is connected
    # before sequence records begin to flow.
    for (group in fifo_by_group) {
        printf "%s", "" > fifo_by_group[group]
        fflush(fifo_by_group[group])
    }

    while ((getline line < mapping_file) > 0) {
        count = split(line, fields, "\t")
        if (count != 2 || !(fields[2] in fifo_by_group)) {
            print "PHAP_ERROR\tinvalid assignment mapping row: " line > "/dev/stderr"
            exit 2
        }
        if (fields[1] in assignment) {
            print "PHAP_ERROR\tduplicate assignment ID: " fields[1] > "/dev/stderr"
            exit 2
        }
        assignment[fields[1]] = fields[2]
        requested++
    }
    close(mapping_file)
}

{
    input_records++
    if (NF != 3) {
        print "PHAP_ERROR\tmalformed seqkit fx2tab record at " input_records > "/dev/stderr"
        parse_errors++
        next
    }

    header = $1
    read_id = header
    sub(/[[:space:]].*$/, "", read_id)
    sub(/^[@>]/, "", read_id)
    sub(/\/[12]$/, "", read_id)
    if (!(read_id in assignment)) {
        if (progress_every > 0 && input_records % progress_every == 0) {
            print "PHAP_PROGRESS\t" input_records > "/dev/stderr"
            fflush("/dev/stderr")
        }
        next
    }
    group = assignment[read_id]
    if (read_id in found) {
        print "PHAP_ERROR\tduplicate assigned FASTQ ID: " read_id > "/dev/stderr"
        duplicate_errors++
        next
    }
    found[read_id] = 1

    filtered_record = (length($2) < min_length)
    if (!filtered_record && min_quality > 0) {
        quality_sum = 0
        quality_length = length($3)
        for (quality_index = 1; quality_index <= quality_length; quality_index++) {
            quality_character = substr($3, quality_index, 1)
            if (!(quality_character in phred_by_character)) {
                print "PHAP_ERROR\tinvalid FASTQ quality character for " read_id > "/dev/stderr"
                parse_errors++
                filtered_record = 1
                break
            }
            quality_sum += phred_by_character[quality_character]
        }
        if (!filtered_record && (quality_length == 0 || quality_sum / quality_length < min_quality)) {
            filtered_record = 1
        }
    }
    if (filtered_record) {
        filtered++
    } else {
        fifo = fifo_by_group[group]
        print "@" header > fifo
        print $2 > fifo
        print "+" > fifo
        print $3 > fifo
        written++
        group_written[group]++
    }

    if (progress_every > 0 && input_records % progress_every == 0) {
        print "PHAP_PROGRESS\t" input_records > "/dev/stderr"
        fflush("/dev/stderr")
    }
}

END {
    for (group in fifo_by_group) {
        close(fifo_by_group[group])
    }

    missing = 0
    missing_examples = ""
    for (read_id in assignment) {
        if (!(read_id in found)) {
            missing++
            if (missing <= 10) {
                missing_examples = missing_examples (missing_examples == "" ? "" : ",") read_id
            }
        }
    }

    print "metric", "value" > summary_file
    print "input_records", input_records + 0 >> summary_file
    print "requested_assigned_reads", requested + 0 >> summary_file
    print "written_reads", written + 0 >> summary_file
    print "filtered_assigned_reads", filtered + 0 >> summary_file
    print "missing_reads", missing + 0 >> summary_file
    print "missing_examples", missing_examples >> summary_file
    print "duplicate_errors", duplicate_errors + 0 >> summary_file
    print "parse_errors", parse_errors + 0 >> summary_file
    for (group in fifo_by_group) {
        print "group_reads:" group, group_written[group] + 0 >> summary_file
    }
    close(summary_file)

    if (duplicate_errors || parse_errors) {
        exit 3
    }
}

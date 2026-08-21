import csv
import pysam


def _infer_tag_value(value):
    """
    Infer a BAM tag value and value_type from a TSV cell string, trying int,
    then float, then falling back to string.
    """
    try:
        return int(value), 'i'
    except ValueError:
        pass
    try:
        return float(value), 'f'
    except ValueError:
        pass
    return value, 'Z'


def _parse_tag_mapping(tsv_path, read_id_column='read_id'):
    """
    Parse a TSV file mapping read_id to tag values. Every column other than
    read_id_column is treated as a 2-letter BAM tag to apply.

    Returns (mapping, tag_columns), where mapping is a dict:
    read_id -> {tag: (value, value_type)}. Blank cells are omitted so that
    tag is left off for that read.
    """
    with open(tsv_path, 'r', newline='') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        columns = reader.fieldnames
        if not columns:
            raise ValueError("The TSV file is empty or improperly formatted.")
        if read_id_column not in columns:
            raise ValueError(f"Column '{read_id_column}' not found in TSV. Available columns: {', '.join(columns)}")

        tag_columns = [c for c in columns if c != read_id_column]
        if not tag_columns:
            raise ValueError(f"TSV must have at least one tag column besides '{read_id_column}'.")
        for tag in tag_columns:
            if len(tag) != 2:
                raise ValueError(f"Invalid tag column '{tag}'. BAM tag columns must be exactly 2 characters.")

        mapping = {}
        for row in reader:
            read_id = row[read_id_column]
            tag_values = {}
            for tag in tag_columns:
                value = row[tag]
                if value is None or value == '':
                    continue
                tag_values[tag] = _infer_tag_value(value)
            mapping[read_id] = tag_values
    return mapping, tag_columns


def _output_mode_from_path(path):
    lower = path.lower()
    if lower.endswith(".cram"):
        return "wc"
    if lower.endswith(".sam"):
        return "w"
    return "wb"


def add_tags_file(input_path, output_path, tsv_path):
    """
    Tag reads in a BAM/CRAM/SAM file with values from a TSV file, keyed by read_id.

    Args:
        input_path: Path to input BAM/CRAM/SAM.
        output_path: Path to output BAM/CRAM/SAM.
        tsv_path: Path to TSV file with a read_id column and one column per
            tag to apply. Each non-read_id column name must be a 2-letter
            BAM tag (e.g., XX, RG).
    """
    mapping, _tags = _parse_tag_mapping(tsv_path)

    with pysam.AlignmentFile(input_path, "rb", check_sq=False) as infile:
        mode = _output_mode_from_path(output_path)
        reference_filename = None
        if mode == "wc":
            reference_filename = getattr(infile, "reference_filename", None)
            if reference_filename is None:
                raise ValueError("CRAM output requires a reference. Use BAM/SAM output or input CRAM with reference.")

        with pysam.AlignmentFile(
            output_path,
            mode,
            template=infile,
            reference_filename=reference_filename,
        ) as outfile:
            total_reads = 0
            tagged_count = 0
            for read in infile:
                total_reads += 1
                # Only primary alignments are tagged, consistent with bamurai's
                # design principle of ignoring secondary/supplementary alignments.
                if not (read.is_secondary or read.is_supplementary):
                    tag_values = mapping.get(read.query_name)
                    if tag_values:
                        for tag, (value, value_type) in tag_values.items():
                            read.set_tag(tag, value, value_type=value_type)
                        tagged_count += 1
                outfile.write(read)

    print(f"Tagged {tagged_count} of {total_reads} reads written to {output_path}")


def add_tags(args):
    """
    CLI wrapper for add_tags_file.
    """
    return add_tags_file(args.bam, args.output, args.tsv)

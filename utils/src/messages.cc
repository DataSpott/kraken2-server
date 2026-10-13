#include "messages.h"

bool SequenceRequestToSequence(
    const kraken2proto::Kraken2SequenceRequest &req, kraken2::Sequence &seq, kraken2::Sequence &seq2)
{
    switch (req.format())
    {
    case kraken2proto::Kraken2SequenceRequest_SequenceFormat_FORMAT_AUTO_DETECT:
        seq.format = kraken2::SequenceFormat::FORMAT_AUTO_DETECT;
        break;
    case kraken2proto::Kraken2SequenceRequest_SequenceFormat_FORMAT_FASTQ:
        seq.format = kraken2::SequenceFormat::FORMAT_FASTQ;
        break;
    case kraken2proto::Kraken2SequenceRequest_SequenceFormat_FORMAT_FASTA:
        seq.format = kraken2::SequenceFormat::FORMAT_FASTA;
        break;
    }

    seq.header = std::string(req.id());
    seq.seq = std::string(req.seq());
    seq.quals = std::string(req.quals());

    if (req.paired()){
        seq2.format = seq.format;
        seq2.header = std::string(req.id2());
        seq2.seq = std::string(req.seq2());
        seq2.quals = std::string(req.quals2());
    }

    return true;
}

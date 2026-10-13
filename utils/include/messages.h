#include "Kraken2.grpc.pb.h"
#include "seqreader.h"

bool SequenceRequestToSequence(
    const kraken2proto::Kraken2SequenceRequest &req, kraken2::Sequence &seq, kraken2::Sequence &seq2);

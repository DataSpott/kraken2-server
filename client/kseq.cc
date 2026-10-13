#include <string>
#include <iostream>
#include <sysexits.h>
#include <cstdlib>
#include "kseq.cc.h"


FastReader::FastReader(std::string filename, const std::string &mates_file)
{
    m_filename = filename;
    m_mates_filename = mates_file;
    m_paired = !mates_file.empty();
    FILE *instream = (filename == "-") ? stdin : fopen(filename.c_str(), "r");
    if (!instream) {
        throw std::runtime_error("Cannot open input: " + filename);
    }

    m_fp = gzdopen(fileno(instream), "r");
    if (!m_fp) {
        if (instream != stdin) {
            fclose(instream);
        }
        throw std::runtime_error("gzdopen failed: " + filename);
    }

    m_seq = kseq_init(m_fp);

    m_fp2 = nullptr;
    m_seq2 = nullptr;
    if (m_paired) {
        FILE *instream2 = fopen(mates_file.c_str(), "r");
        if (!instream2) {
            kseq_destroy(m_seq);
            gzclose(m_fp);
            std::cerr << "Cannot open mates input: " + mates_file << std::endl;
            exit(EX_DATAERR);
        }
        m_fp2 = gzdopen(fileno(instream2), "r");
        if (!m_fp2) {
            fclose(instream2);
            kseq_destroy(m_seq);
            gzclose(m_fp);
            std::cerr << "gzdopen failed: " + mates_file << std::endl;
            exit(EX_DATAERR);
        }
        m_seq2 = kseq_init(m_fp2);
    }
}


FastReader::~FastReader()
{
    kseq_destroy(m_seq);
    gzclose(m_fp);
    if (m_paired) {
        kseq_destroy(m_seq2);
        gzclose(m_fp2);
    }
}


int FastReader::read(Kraken2SequenceRequest& rec) {
    int rtn;
    if ((rtn = kseq_read(m_seq)) < 0) {
        if (m_paired && kseq_read(m_seq2) >= 0) {
            std::cerr << "Mates file longer than sequence file (" + m_mates_filename + ")" << std::endl;
            exit(EX_DATAERR);
        }
        rec.Clear();
        return rtn;
    }
    else {
        std::string header;
        m_seq->qual.l == 0 ? header.append(">") : header.append("@");
        header.append(m_seq->name.s);
        if(m_seq->comment.l > 0) {
            header.append(" ");
            header.append(m_seq->comment.s);
        }

        std::string str_rep;
        str_rep.append(header);
        str_rep.append("\n");

        rec.set_id(m_seq->name.s);
        rec.set_seq(m_seq->seq.s);
        rec.set_header(header);
        if (m_seq->qual.l == 0)
        {
            rec.set_format(Kraken2SequenceRequest::FORMAT_FASTA);
            rec.set_quals("");
            str_rep.append(rec.seq());
            str_rep.append("\n");
        }
        else
        {
            rec.set_format(Kraken2SequenceRequest::FORMAT_FASTQ);
            rec.set_quals(m_seq->qual.s);
            str_rep.append(rec.seq());
            str_rep.append("\n+\n");
            str_rep.append(rec.quals());
            str_rep.append("\n");
        }
        rec.set_str_representation(str_rep);

        if (m_paired) {
            if (kseq_read(m_seq2) < 0) {
                rec.Clear();
                std::cerr << "Mates file ended before sequence file (" + m_mates_filename + ")" << std::endl;
                exit(EX_DATAERR);
            }
            std::string name1(m_seq->name.s);
            std::string name2(m_seq2->name.s);
            auto strip_mate_suffix = [](const std::string &n) {
                if (n.size() > 2 && n[n.size() -2] == '/' && (n.back() == '1' || n.back() == '2'))
                    return n.substr(0, n.size() -2);
                return n;
            };
            if (strip_mate_suffix(name1) != strip_mate_suffix(name2)) {
                rec.Clear();
                std::cerr << "Unpaired mates: '" + name1 + "' vs '" + name2 + "'" << std::endl;
                exit(EX_DATAERR);
            }
            std::string header2;
            m_seq2->qual.l == 0 ? header2.append(">") : header2.append("@");
            header2.append(m_seq2->name.s);
            if (m_seq2->comment.l > 0) {
                header2.append(" ");
                header2.append(m_seq2->comment.s);
            }
            rec.set_paired(true);
            rec.set_id2(m_seq2->name.s);
            rec.set_header2(header2);
            rec.set_seq2(m_seq2->seq.s);
            if (m_seq2->qual.l == 0) {
                rec.set_quals2("");
            } else {
                rec.set_quals2(m_seq2->qual.s);
            }
        }
    }
    return rtn;
}


int FastReader::read(std::vector<Kraken2SequenceRequest> &seqs, int batch_size)
{
    int rtn = 0;
    seqs.reserve(batch_size);
    for(size_t i=0; i<batch_size; ++i) {
        Kraken2SequenceRequest rec;
        if (read(rec) >= 0) {
            rtn++;
            seqs.push_back(std::move(rec));
        }
        else {
            break;
        }
    }
    return rtn;
}


int FastReader::read_all(std::vector<Kraken2SequenceRequest> &seqs)
{
    int rtn = 0;
    while (true) {
        Kraken2SequenceRequest rec;
        if (read(rec) >= 0) {
            rtn++;
            seqs.push_back(std::move(rec));
        }
        else {
            break;
        }
    }
    return rtn;
}

/**
 * MALVA - genotyping by Mapping-free ALternate-allele detection of known VAriants
 * Copyright (C) 2019  Giulia Bernardini, Luca Denti, Marco Previtali
 *
 * This file is part of MALVA.
 *
 * MALVA is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MALVA is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MALVA; see the file LICENSE. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <utility>

#include <unordered_map>
#include <set>
#include <string>
#include <vector>

#include <math.h>
#include <zlib.h>

#include "hts_log.h"
#include "kmc_file.h"
#include "kseq.h"
#include "vcf.h"

#include "argument_parser.hpp"
#include "bloom_filter.hpp"
#include "var_block.hpp"
#include "kmap.hpp"

using namespace std;

auto start_t = chrono::high_resolution_clock::now();

void pelapsed(const string &s = "", const bool rollback = false) {
  auto now_t = chrono::high_resolution_clock::now();
  cerr << "[malva-geno/" << s << "] Time elapsed "
	    << chrono::duration_cast<chrono::milliseconds>(now_t - start_t).count()/1000 << "s";
  if(rollback) cerr << "\r";
  else cerr << endl;
}

KSEQ_INIT(gzFile, gzread)

/**
 * Method to add kmers to the bloom filter
 **/
void add_kmers_to_bf(BF &bf, KMAP &ref_bf, const VK_GROUP &kmers) {
  for (const auto &v : kmers) {
    // For each variant
    for (const auto &p : v.second) {
      // For each allele of the variant
      for (const auto &Ks : p.second) {
        // For each list of kmers of the allele
        for (const auto &kmer : Ks) {
          // For each kmer in the kmer list
	  if(p.first == 0)
            ref_bf.add_key(kmer.c_str());
	  else
            bf.add_key(kmer.c_str());
        }
      }
    }
  }
}

/**
 * Method to compute and store the coverages of the alleles of the
 * variants of a var_block. It uses the coverages stored in the bloom
 * filters/map.
 **/
void set_coverages(BF &bf, KMAP &ref_bf, VB &vb, const VK_GROUP &kmers/*, const float &cap*/) {
  for (const auto &var : kmers) {
    // For each variant
    Variant v = vb.get_variant(var.first);
    for (const auto &p : var.second) {
      uint allele_cov = 0;
      for (const auto &Ks : p.second) {
        uint curr_cov = 0;
        int n = 0; // Number of kmers in the signature
        for (const auto &kmer : Ks) {
          int w = 0;
          if(p.first == 0)
            w = ref_bf.get_count(kmer.c_str());
          else
            w = bf.get_count(kmer.c_str());
          if(w>0) { //maybe useless
            curr_cov = (curr_cov * n + w) / (n + 1);
            ++n;
          }
        }
        if(curr_cov > allele_cov)
          allele_cov = curr_cov;
      }
      // we can now set the allele coverage
      vb.set_variant_coverage(var.first, p.first, allele_cov);
    }
  }
}

/**
 * Method to clean and print VCF header. It adds GT and GQ FORMAT,
 * removes all samples, and adds donor sample.
 **/
void print_cleaned_header(bcf_hdr_t *vcf_header, const bool verbose) {
  // Adding format fields - if already present, they won't be added
  bcf_hdr_append(vcf_header, "##FORMAT=<ID=GT,Number=1,Type=String,"
                 "Description=\"Genotype\">");
  bcf_hdr_append(vcf_header, "##FORMAT=<ID=GQ,Number=1,Type=Integer,"
                 "Description=\"Genotype Quality\">");

  if(verbose) {
    bcf_hdr_append(vcf_header, "##INFO=<ID=COVS,Number=R,Type=Integer,"
		   "Description=\"Allele coverages\">");
    bcf_hdr_append(vcf_header, "##INFO=<ID=GTS,Number=.,Type=String,"
		   "Description=\"Genotypes Likelihood\">");
  }

  // Adding donor sample and removing all other samples
  const char *new_sample = "DONOR";
  bcf_hdr_add_sample(vcf_header, new_sample);
  bcf_hdr_sync(vcf_header);
  bcf_hdr_set_samples(vcf_header, new_sample, 0);

  // Formatting and printing header
  kstring_t htxt = {0, 0, 0};
  bcf_hdr_format(vcf_header, 0, &htxt);
  cout << htxt.s;
  free(htxt.s);
}

// ---------------------------------------------------------------------------

int main(int argc, char *argv[]) {
  hts_set_log_level(HTS_LOG_OFF);

  parse_arguments(argc, argv);

  // STEP 0: open and check input files
  gzFile fasta_in = gzopen(opt::fasta_path.c_str(), "r");
  kseq_t *reference = kseq_init(fasta_in);

  htsFile *vcf = bcf_open(opt::vcf_path.c_str(), "r");
  bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);
  int is_file_flag = 0;
  if(opt::samples != "-")
    is_file_flag = 1;
  int set_samples_code = bcf_hdr_set_samples(vcf_header, opt::samples.c_str(), is_file_flag);
  if(set_samples_code != 0) {
    cerr << "ERROR: VCF samples subset (code: " << set_samples_code << ")" << endl;
    return 1;
  }
  bcf1_t *vcf_record = bcf_init();

  CKMCFile kmer_db;
  if (!kmer_db.OpenForListing(opt::kmc_sample_path)) {
    cerr << "ERROR: cannot open " << opt::kmc_sample_path << endl;
    return 1;
  }

  // References are stored in a map
  pelapsed("Reference parsing");
  unordered_map<string, string> refs;
  int l;
  while ((l = kseq_read(reference)) >= 0) {
    string id = reference->name.s;
    if(id.compare(0,3,"chr") == 0 && opt::strip_chr) {
      id = id.substr(3);
    }
    string seq (reference->seq.s);
    refs[id] = seq;
  }
  pelapsed("Reference processed");

  // STEP 1: add VCF kmers to bloom filter
  pelapsed("VCF parsing (Bloom Filter construction)");
  BF bf(opt::bf_size);
  KMAP ref_bf;
  // BF context_bf(opt::bf_size);

  vector<string> used_seq_names;
  VB vb(opt::k, opt::error_rate);
  string last_seq_name = "";

  int i = 0;
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_STR);
    Variant v(vcf_header, vcf_record, opt::freq_key, opt::uniform);
    ++i;
    if(i%5000 == 0) {
      string log_line = "Processed " + to_string(i) + " variants";
      pelapsed(log_line);
    }

    // In the first iteration, we set last_seq_name
    if(last_seq_name.size() == 0) {
      last_seq_name = v.seq_name;
      used_seq_names.push_back(last_seq_name);
    }

    // We do not consider variants with <CN> or not present in
    // considered samples, i.e. 0|0 for all samples
    if (!v.has_alts or !v.is_present)
      continue;

    if (vb.empty()) {
      vb.add_variant(v);
      continue;
    }

    if (!vb.is_near_to_last(v) || last_seq_name != v.seq_name) {
      /***
       * 1. extract k-mers
       * 2. add k-mers to BF
       * 3. clear block
       * 4. set new reference
       ***/
      VK_GROUP kmers = vb.extract_kmers(refs[last_seq_name], opt::haploid);
      add_kmers_to_bf(bf, ref_bf, kmers);
      vb.clear();
      if(last_seq_name != v.seq_name) {
        last_seq_name = v.seq_name;
        used_seq_names.push_back(last_seq_name);
      }
    }
    vb.add_variant(v);
  }
  if (!vb.empty()) {
    /***
     * 1. extract k-mers
     * 2. add k-mers to BF
     * 3. clear block
     ***/
    VK_GROUP kmers = vb.extract_kmers(refs[last_seq_name], opt::haploid);
    add_kmers_to_bf(bf, ref_bf, kmers);
    vb.clear();
  }
  string log_line = "Processed " + to_string(i) + " variants";
  pelapsed(log_line);

  bcf_hdr_destroy(vcf_header);
  bcf_destroy(vcf_record);
  bcf_close(vcf);

  bf.switch_mode();

  pelapsed("BF creation complete");

  // pelapsed("Reference BF construction");
  // for(const auto &seq_name : used_seq_names) {
  //   string reference = refs[seq_name];
  //   string ref_ksub(reference, (opt::ref_k - opt::k) / 2, opt::k);
  //   string context(reference, 0, opt::ref_k);
  //   transform(ref_ksub.begin(), ref_ksub.end(), ref_ksub.begin(), ::toupper);
  //   transform(context.begin(), context.end(), context.begin(), ::toupper);
  //   if (bf.test_key(ref_ksub.c_str()))
  //     context_bf.add_key(context.c_str());
  //   for (uint p = opt::ref_k; p < reference.size(); ++p) {
  //     char c1 = toupper(reference[p]);
  //     context.erase(0, 1);
  //     context += c1;
  //     char c2 = toupper(reference[p - (opt::ref_k - opt::k) / 2]);
  //     ref_ksub.erase(0, 1);
  //     ref_ksub += c2;
  //     if (bf.test_key(ref_ksub.c_str()))
  //       context_bf.add_key(context.c_str());
  //   }
  // }
  // pelapsed("Reference BF creation complete");

  // context_bf.switch_mode();

  // STEP 2: test variants present in read sample
  pelapsed("KMC output processing");
  uint32 klen, mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64 tot_kmers, max_c;
  kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c, tot_kmers);
  CKmerAPI kmer_obj(klen);

  char kmer[opt::k + 1];
  while (kmer_db.ReadNextKmer(kmer_obj, counter)) {
    kmer_obj.to_string(kmer);
    transform(kmer, kmer + opt::k, kmer, ::toupper);
    // char kmer[opt::k + 1];
    // strncpy(kmer, context + ((opt::ref_k - opt::k) / 2), opt::k);
    // kmer[opt::k] = '\0';
    ref_bf.increment(kmer, counter);
    // if (!context_bf.test_key(context)) {
    bf.increment(kmer, counter);
    // }
  }

  pelapsed("BF weights created");

  // STEP 3: check if variants in vcf are covered enough
  vcf = bcf_open(opt::vcf_path.c_str(), "r");
  vcf_header = bcf_hdr_read(vcf);
  set_samples_code = bcf_hdr_set_samples(vcf_header, NULL, 0);
  print_cleaned_header(vcf_header, opt::verbose);
  bcf_hdr_destroy(vcf_header);
  bcf_close(vcf);

  vcf = bcf_open(opt::vcf_path.c_str(), "r");
  vcf_header = bcf_hdr_read(vcf);
  set_samples_code = bcf_hdr_set_samples(vcf_header, opt::samples.c_str(), is_file_flag);
  vcf_record = bcf_init();

  pelapsed("VCF parsing and genotyping");
  i = 0;
  last_seq_name = "";
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_STR);
    Variant v(vcf_header, vcf_record, opt::freq_key, opt::uniform);
    ++i;
    if(i%1000 == 0) {
      string log_line = "Processed " + to_string(i) + " variants";
      pelapsed(log_line);
    }
    // In the first iteration, we set last_seq_name
    if(last_seq_name.size() == 0)
      last_seq_name = v.seq_name;

    // In this step, we must consider the variants not present in the
    // samples: their genotype is 0/0
    if (!v.has_alts)
      continue;

    if (vb.empty()) {
      vb.add_variant(v);
      continue;
    }

    if (!vb.is_near_to_last(v) || last_seq_name != v.seq_name) {
      /***
       * 1. extract k-mers
       * 2. check if variants are covered
       * 3. output covered variants
       * 4. clear block
       * 5. set new reference
       ***/
      VK_GROUP kmers = vb.extract_kmers(refs[last_seq_name], opt::haploid);
      set_coverages(bf, ref_bf, vb, kmers);
      vb.genotype(opt::max_coverage, opt::haploid);
      vb.output_variants(opt::haploid, opt::verbose);
      cout.flush();
      vb.clear();
      if(last_seq_name != v.seq_name)
        last_seq_name = v.seq_name;
    }
    vb.add_variant(v);
  }
  if (!vb.empty()) {
    /***
     * 1. extract k-mers
     * 2. check if variants are covered
     * 3. output covered variants
     * 4. clear block
     ***/
    VK_GROUP kmers = vb.extract_kmers(refs[last_seq_name], opt::haploid);
    set_coverages(bf, ref_bf, vb, kmers);
    vb.genotype(opt::max_coverage, opt::haploid);
    vb.output_variants(opt::haploid, opt::verbose);
    cout.flush();
    vb.clear();
  }
  log_line = "Processed " + to_string(i) + " variants";
  pelapsed(log_line);

  bcf_hdr_destroy(vcf_header);
  bcf_destroy(vcf_record);
  bcf_close(vcf);

  kseq_destroy(reference);
  gzclose(fasta_in);

  cout.flush();

  pelapsed("Execution completed");

  return 0;
}

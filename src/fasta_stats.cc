#include "open_compressed.h"	/* close_compressed(), open_compressed(), pfgets() */
#include "pretty_print.h"	/* pretty_print() */
#include <errno.h>	/* errno */
#include <list>		/* list<> */
#include <stdio.h>	/* EOF, fprintf(), printf(), stderr */
#include <stdlib.h>	/* atol(), getopt(), optarg, optind */
#include <string.h>	/* strerror() */
#include <string>	/* string */
#include <unistd.h>	/* exit() */
#include <utility>	/* pair<> */
#include <vector>	/* vector<> */

static long opt_n_cutoff, opt_N_cutoff, opt_big_size, opt_min_length;
static long opt_min_gap_length;

class Contig {
    public:
	long total_size;	// includes n's
	long n;
	Contig(void) : total_size(0), n(0) { }
	~Contig(void) { }
	bool operator<(const Contig &__a) const {
		return size() > __a.size();
	}
	long size(void) const {
		return total_size - n;
	}
	bool empty(void) {
		if (total_size == n) {
			total_size = n = 0;
			return 1;
		}
		return 0;
	}
};

class Scaffold {
    public:
	long total_size, contig_size, gap_size, captured_gap_size;
	std::list<Contig> contigs;
	std::list<long> gaps;
	Scaffold(void) : total_size(0), contig_size(0), gap_size(0), captured_gap_size(0), contigs(1, Contig()) { }
	~Scaffold(void) { }
	bool operator<(const Scaffold &__a) const {
		return total_size > __a.total_size;
	}
	void calc_contig_sizes(void) {
		gap_size = total_size;
		std::list<Contig>::const_iterator __b = contigs.begin();
		std::list<Contig>::const_iterator __end_b = contigs.end();
		for (; __b != __end_b; ++__b) {
			contig_size += __b->size();
			gap_size -= __b->total_size;
		}
		std::list<long>::const_iterator __c = gaps.begin();
		std::list<long>::const_iterator __end_c = gaps.end();
		for (; __c != __end_c; ++__c) {
			captured_gap_size += *__c;
		}
	}
	bool empty(void) {
		if (contigs.back().empty() && contigs.size() > 1) {
			contigs.pop_back();
			gaps.pop_back();
		}
		if (contigs.size() == 1 && contigs.back().total_size == 0) {
			total_size = 0;
			return 1;
		}
		return 0;
	}
};

class Stats {
    public:
	long scaffolds, contigs, s_size, c_size;
	Stats(void) : scaffolds(0), contigs(0), s_size(0), c_size(0) { }
	~Stats(void) { }
	void add(const Scaffold &__a) {
		++scaffolds;
		contigs += __a.contigs.size();
		s_size += __a.total_size;
		c_size += __a.contig_size;
	}
	void add(const Stats &__a) {
		scaffolds += __a.scaffolds;
		contigs += __a.contigs;
		s_size += __a.s_size;
		c_size += __a.c_size;
	}
	double coverage(void) const {
		return s_size == 0 ? 0 : (double)100 * c_size / s_size;
	}
};

// read in fasta file
static int read_file(char *filename, std::list<Scaffold> &scaffolds) {
	int fd = open_compressed(filename);
	if (fd == -1) {
		fprintf(stderr, "Error: open_compressed %s: %s\n", filename, strerror(errno));
		return 0;
	}
	long n = 0;
	scaffolds.push_back(Scaffold());
	Contig *current = &scaffolds.back().contigs.back();
	std::string line;
	while (pfgets(fd, line) != -1) {
		if (line[0] == '>') {
			if (!scaffolds.back().empty()) {
				scaffolds.push_back(Scaffold());
			}
			current = &scaffolds.back().contigs.back();
			n = 0;
			continue;
		}
		size_t i = line.find_first_of("Nn");
		if (i != 0 && n != 0) {		// handle line carryover
			if (n == opt_N_cutoff) {
				if (!scaffolds.back().empty()) {
					scaffolds.back().total_size -= n;
					scaffolds.push_back(Scaffold());
				}
				current = &scaffolds.back().contigs.back();
			} else if (n >= opt_n_cutoff) {
				if (!current->empty()) {
					scaffolds.back().contigs.push_back(Contig());
					scaffolds.back().gaps.push_back(n);
					current = &scaffolds.back().contigs.back();
				}
			} else {
				current->n += n;
				current->total_size += n;
			}
			n = 0;
		}
		scaffolds.back().total_size += line.size();
		size_t j = 0;
		for (;;) {
			if (i == line.npos) {
				current->total_size += line.size() - j;
				break;
			}
			current->total_size += i - j;
			j = line.find_first_not_of("Nn", i + 1);
			if (j == line.npos) {
				n += line.size() - i;
				break;
			}
			n += j - i;
			if (n == opt_N_cutoff) {
				if (!scaffolds.back().empty()) {
					scaffolds.back().total_size -= n + line.size() - j;
					scaffolds.push_back(Scaffold());
				}
				current = &scaffolds.back().contigs.back();
				scaffolds.back().total_size = line.size() - j;
			} else if (n >= opt_n_cutoff) {
				if (!current->empty()) {
					scaffolds.back().contigs.push_back(Contig());
					scaffolds.back().gaps.push_back(n);
					current = &scaffolds.back().contigs.back();
				}
			} else {
				current->total_size += n;
				current->n += n;
			}
			n = 0;
			i = line.find_first_of("Nn", j + 1);
		}
	}
	close_compressed(fd);
	if (scaffolds.back().empty()) {
		scaffolds.pop_back();
	}
	std::list<Scaffold>::iterator a = scaffolds.begin();
	std::list<Scaffold>::iterator end_a = scaffolds.end();
	for (; a != end_a; ++a) {
		a->calc_contig_sizes();
	}
	scaffolds.sort();
	return 1;
}

static void scaled_print(long i) {
	if (i > 999999) {
		printf("%.1f MB", (double)i / 1000000);
	} else if (i > 999) {
		printf("%.1f KB", (double)i / 1000);
	} else {
		printf("%ld", i);
	}
}

static void print_stats(const std::list<Scaffold> &scaffolds) {
	// make contig list
	std::list<Contig> contigs;
	std::list<Scaffold>::const_iterator a = scaffolds.begin();
	std::list<Scaffold>::const_iterator end_a = scaffolds.end();
	for (; a != end_a; ++a) {
		contigs.insert(contigs.end(), a->contigs.begin(), a->contigs.end());
	}
	contigs.sort();
	// collect statistics
	long scaffold_seq = 0;
	long contig_seq = 0;
	long gap_seq = 0;
	long big_scaffolds = 0;
	long big_scaffold_seq = 0;
	for (a = scaffolds.begin(); a != end_a; ++a) {
		scaffold_seq += a->total_size;
		contig_seq += a->contig_size;
		gap_seq += a->gap_size;
		if (a->total_size > opt_big_size) {
			++big_scaffolds;
			big_scaffold_seq += a->total_size;
		}
	}
	long n50_scafs = 0;
	long x = 0;
	for (a = scaffolds.begin(); a != end_a && 2 * x < scaffold_seq; ++a) {
		++n50_scafs;
		x += a->total_size;
	}
	if (a != scaffolds.begin()) {
		--a;
	}
	long n50_contigs = 0;
	std::list<Contig>::const_iterator b = contigs.begin();
	std::list<Contig>::const_iterator end_b = contigs.end();
	for (x = 0; b != end_b && 2 * x < contig_seq; ++b) {
		++n50_contigs;
		x += b->size();
	}
	if (b != contigs.begin()) {
		--b;
	}
	// print statistics
	printf("Main genome scaffold total: %lu\n", scaffolds.size());
	printf("Main genome contig total:   %lu\n", contigs.size());
	printf("Main genome scaffold sequence total: ");
	scaled_print(scaffold_seq);
	printf("\nMain genome contig sequence total:   ");
	scaled_print(contig_seq);
	printf(" (-> %4.1f%% gap)\n", (double)100 * gap_seq / scaffold_seq);
	printf("Main genome scaffold N/L50: %ld/", n50_scafs);
	scaled_print(a->total_size);
	printf("\nMain genome contig N/L50:   %ld/", n50_contigs);
	scaled_print(b->size());
	printf("\nNumber of scaffolds > %ld KB: %ld\n", opt_big_size / 1000, big_scaffolds);
	printf("%% main genome in scaffolds > %ld KB: %4.1f%%\n", opt_big_size / 1000, (double)100 * big_scaffold_seq / scaffold_seq);
}

static void print_size(const int lengths[], int i, int one_less = 0, int right_justify = 0) {
	int x = lengths[i] - (one_less ? 1 : 0);
	if ((i % 3) == 2) {	// 250... line
		if (x >= 1000000) {
			printf("%3.1f mb", (double)x / 1000000);
		} else if (x >= 10000) {
			printf("%3d kb", x / 1000);
		} else if (x >= 1000) {
			printf("%3.1f kb", (double)x / 1000);
		} else if (right_justify) {
			printf("   %3d", x);
		} else {
			printf("%3d   ", x);
		}
	} else {
		if (x >= 1000000) {
			printf("%3d mb", x / 1000000);
		} else if (x >= 1000) {
			printf("%3d kb", x / 1000);
		} else if (right_justify) {
			printf("   %3d", x);
		} else {
			printf("%3d   ", x);
		}
	}
}

static void print_coverage(const std::list<Scaffold> &scaffolds) {
	const int lengths[] = {       0,
				    100,     250,     500,
				   1000,    2500,    5000,
				  10000,   25000,   50000,
				 100000,  250000,  500000,
				1000000, 2500000, 5000000
	};
	const size_t lengths_size = sizeof(lengths) / sizeof(int);
	size_t i;
	std::vector<Stats> stats(lengths_size);
	std::list<Scaffold>::const_iterator a = scaffolds.begin();
	std::list<Scaffold>::const_iterator end_a = scaffolds.end();
	for (; a != end_a; ++a) {
		for (i = lengths_size - 1; a->total_size < lengths[i]; --i) { }
		stats[i].add(*a);
	}
	for (i = lengths_size - 1; i != 0; --i) {
		stats[i - 1].add(stats[i]);
	}
	printf(" Minimum    Number    Number     Total        Total     Scaffold\n");
	printf("Scaffold      of        of      Scaffold      Contig     Contig\n");
	printf(" Length   Scaffolds  Contigs     Length       Length    Coverage\n");
	printf("--------  ---------  -------  -----------  -----------  --------\n");
	int found_end = 0;
	for (i = 0; i != lengths_size; ++i) {
		const int &x = lengths[i];
		const Stats &b = stats[i];
		if (x == 0) {
			printf("    All");
		} else if (x < opt_min_length || (found_end && b.scaffolds == 0)) {
			continue;
		} else {
			printf(" ");
			print_size(lengths, i);
		}
		if (b.scaffolds == 0) {
			found_end = 1;
		}
		printf("   %7s    %7s  %11s  %11s   %6.2f%%\n", pretty_print(b.scaffolds).c_str(), pretty_print(b.contigs).c_str(), pretty_print(b.s_size).c_str(), pretty_print(b.c_size).c_str(), b.coverage());
	}
}

static void print_captured_gaps(const std::list<Scaffold> &scaffolds) {
	const int lengths[] = {       1,
				    100,     200,     300,
				    400,     500,     600,
				    700,     800,     900,
				   1000,    2500,    5000,
				  10000,   25000,   50000,
				 100000,  250000,  500000,
				1000000, 2500000, 5000000
	};
	const size_t lengths_size = sizeof(lengths) / sizeof(int);
	size_t i;
	std::vector<std::pair<long, long> > gaps(lengths_size, std::pair<long, long>(0, 0));
	std::list<Scaffold>::const_iterator a = scaffolds.begin();
	std::list<Scaffold>::const_iterator end_a = scaffolds.end();
	for (; a != end_a; ++a) {
		std::list<long>::const_iterator b = a->gaps.begin();
		std::list<long>::const_iterator end_b = a->gaps.end();
		for (; b != end_b; ++b) {
			for (i = lengths_size - 1; *b < lengths[i]; --i) { }
			++gaps[i].first;
			gaps[i].second += *b;
		}
	}
	size_t last_found = 0;
	std::pair<long, long> total_gaps(0, 0);
	std::vector<std::pair<long, long> >::const_iterator b = gaps.begin();
	std::vector<std::pair<long, long> >::const_iterator end_b = gaps.end();
	for (i = 0; b != end_b; ++b, ++i) {
		total_gaps.first += b->first;
		total_gaps.second += b->second;
		if (b->first != 0) {
			last_found = i;
		}
	}
	printf("                  Number                       Percent     Total     Percent\n");
	printf("     Gap            of        Gap      Total    Total       Gap       Total\n");
	printf("  Size Range       Gaps     Lengths     Gaps     Gaps     Lengths    Lengths\n");
	printf("---------------  -------  ----------  -------  -------  -----------  -------\n");
	std::pair<long, long> gap_count(0, 0);
	for (i = 0, ++last_found; i != lengths_size && lengths[i] <= opt_min_gap_length; ++i) { }
	int first_print = i > 1;
	if (i == 1) {
		i = 0;
	} else if (i > 1) {
		i -= 2;
	}
	for (; i != last_found; ++i) {
		const std::pair<long, long> &c = gaps[i];
		gap_count.first += c.first;
		gap_count.second += c.second;
		if (first_print) {
			first_print = 0;
			printf("       < ");
			print_size(lengths, i + 1);
		} else {
			print_size(lengths, i, 0, 1);
			printf(" - ");
			print_size(lengths, i + 1, 1);
		}
		printf("  %7s  %10s  %7s  %6.2f%%  %11s  %6.2f%%\n",
			pretty_print(c.first).c_str(),
			pretty_print(c.second).c_str(),
			pretty_print(gap_count.first).c_str(),
			total_gaps.first == 0 ? 0 : (double)100 * gap_count.first / total_gaps.first,
			pretty_print(gap_count.second).c_str(),
			total_gaps.second == 0 ? 0 : (double)100 * gap_count.second / total_gaps.second
);
	}
}

static void print_usage() {
	fprintf(stderr, "usage: fasta_stats [options] <fasta_file>\n");
	fprintf(stderr, "\t-b ##\tsize of big scaffold cutoff in kb [50]\n");
	fprintf(stderr, "\t-g ##\tminimum size bin to display for gaps [none]\n");
	fprintf(stderr, "\t-l ##\tminimum size bin to display [1000]\n");
	fprintf(stderr, "\t-N ##\texact length of N's to signal scaffold end [-1]\n");
	fprintf(stderr, "\t-n ##\tminimum length of N's to signal contig end [9]\n");
	exit(0);
}

int main(int argc, char **argv) {
	opt_n_cutoff = 9;
	opt_N_cutoff = -1;
	opt_big_size = 50;
	opt_min_length = 1000;
	opt_min_gap_length = -1;
	int c;
	while ((c = getopt(argc, argv, "b:g:l:N:n:")) != EOF) {
		switch (c) {
		    case 'b':
			opt_big_size = atol(optarg);
			break;
		    case 'g':
			opt_min_gap_length = atol(optarg);
			if (opt_min_gap_length < 0) {
				print_usage();
			}
			break;
		    case 'l':
			opt_min_length = atol(optarg);
			if (opt_min_length < 0) {
				print_usage();
			}
			break;
		    case 'N':
			opt_N_cutoff = atol(optarg);
			break;
		    case 'n':
			opt_n_cutoff = atol(optarg);
			if (opt_n_cutoff < 0) {
				print_usage();
			}
			break;
		    default:
			print_usage();
		}
	}
	opt_big_size *= 1000;
	if (optind + 1 != argc) {
		print_usage();
	}
	std::list<Scaffold> scaffolds;
	if (!read_file(argv[optind], scaffolds)) {
		return 1;
	}
	print_stats(scaffolds);
	printf("\n");
	print_coverage(scaffolds);
	if (opt_min_gap_length != -1) {
		printf("\n");
		print_captured_gaps(scaffolds);
	}
	return 0;
}

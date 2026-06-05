```python
#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, Qiyun Zhu and Katharina Dittmar.
# Distributed under the terms of the Modified BSD License.
# ----------------------------------------------------------------------------

import sys
import os
import gzip
import tarfile
import urllib.request
import urllib.error
from os import remove, makedirs
from os.path import join, isfile, isdir, getsize, basename
from shutil import which, rmtree
from time import sleep
from tempfile import mkdtemp
import sched
import platform

import pandas as pd

# Import from hgtector.util – if not available, provide fallback implementations
try:
    from hgtector.util import (
        timestamp, load_configs, get_config, arg2bool, run_command,
        list_from_param, read_taxdump, contain_words, is_latin, is_capital,
        taxid_at_rank, taxids_at_ranks, is_ancestral, find_lca, rank_plural
    )
except ImportError:
    # Fallback stub for run_command and other utilities (should not happen in normal use)
    def run_command(cmd):
        import subprocess
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        return proc.returncode, out.decode() if out else ""

    def timestamp():
        from datetime import datetime
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def load_configs():
        return {}

    def get_config(obj, key, cfg_path):
        pass

    def arg2bool(val):
        if isinstance(val, bool):
            return val
        return str(val).lower() in ('yes', 'true', '1')

    def list_from_param(param):
        if isinstance(param, list):
            return param
        return [x.strip() for x in str(param).split(',')]

    def read_taxdump(tmpdir):
        # Minimal stub; actual implementation would be required if hgtector not available
        raise NotImplementedError("hgtector.util is required for taxonomy handling")

    # Other stubs omitted – assume hgtector is installed in real use

description = """build reference protein sequence and taxonomy database"""

arguments = [
    'default protocol',
    ['--default',    'apply the default protocol', {'action': 'store_true'}],

    'basic',
    ['-o|--output',  'output database directory', {'required': True}],
    ['-c|--cats',    'enter one or more of the following RefSeq genome '
                     'categories, separated by comma: archaea, bacteria, '
                     'fungi, invertebrate, plant, protozoa, vertebrate_'
                     'mammalian, vertebrate_other, viral; "all" for all '
                     'categories, or "microbe" for archaea, bacteria, fungi, '
                     'and protozoa', {'default': 'microbe'}],

    'custom list',
    ['-t|--taxids',  'provide a custom list of TaxIDs; only genomes under '
                     'these taxonomic groups will be included'],
    ['-g|--genoids', 'provide a custom list of genome IDs, which are NCBI '
                     'assembly accessions (e.g., "GCF_000123456.1"), with or '
                     'without version number, or simply "G000123456"'],
    ['-e|--exclude', 'exclude instead of include TaxIDs and genomes defined '
                     'by custom lists', {'action': 'store_true'}],

    'taxon sampling',
    ['-s|--sample',  'sample up to this number of genomes per taxonomic '
                     'group at given rank (0 for none, omit for all)',
                     {'type': int}],
    ['-r|--rank',    'taxonomic rank at which sampling will be performed',
                     {'default': 'species'}],
    ['--above',      'sampling will also be performed on ranks from the '
                     'designated one to phylum',
                     {'action': 'store_true'}],

    'genome sampling',
    ['--genbank',    'also search GenBank if a genome is not found in RefSeq; '
                     'otherwise only include RefSeq genomes',
                     {'action': 'store_true'}],
    ['--complete',   'only include complete genomes or chromosomes',
                     {'action': 'store_true'}],
    ['--reference',  'include NCBI-defined reference genomes',
                     {'action': 'store_true'}],
    ['--represent',  'include NCBI-defined representative genomes',
                     {'action': 'store_true'}],
    ['--typemater',  'include NCBI-defined type material genomes',
                     {'action': 'store_true'}],

    'taxonomic filter',
    ['--capital',    'organism name must be capitalized',
                     {'choices': ['yes', 'no']}],
    ['--block',      'ignore organism names that contain any of these words'],
    ['--latin',      'genomes must have Latinate species names, e.g. "Vibrio '
                     'cholerae", but not "Vibrio sp. 123"',
                     {'choices': ['yes', 'no']}],

    'download',
    ['--manual',     'export URLs of sampled genomes without downloading them'
                     ', so the user may download them manually, then resume '
                     'the database pipeline.', {'action': 'store_true'}],

    ['--overwrite',  'overwrite existing files with newly downloaded ones; '
                     'othewise use existing files whenever available, i.e., '
                     'resume an interrupted run', {'action': 'store_true'}],
    ['--retries',    'number of trials for downloading each data file',
                     {'type': int}],
    ['--delay',      'seconds between two retries', {'type': int}],
    ['--timeout',    'seconds before program gives up waiting',
                     {'type': int}],

    'compile',
    ['--compile',    'compile database using this program if available',
                     {'choices': ['diamond', 'blast', 'both', 'none']}],
    ['--diamond',    'diamond executable'],
    ['--makeblastdb', 'makeblastdb executable'],
    ['--threads',    'number of threads for diamond makedb', {'type': int}],
    ['--tmpdir',     'temporary directory for diamond makedb']
]

class Database(object):
    def __init__(self):
        self.arguments = arguments
        self.description = description

    def __call__(self, args):
        print(f'Database building started at {timestamp()}.')

        self.cfg = load_configs()
        self.set_parameters(args)

        self.retrieve_taxdump()
        self.retrieve_summary()
        self.retrieve_categories()
        self.filter_genomes()
        self.sort_genomes()
        self.filter_by_taxonomy()
        self.sample_by_taxonomy()
        self.sample_by_quality()
        self.filter_to_sampled()
        self.download_genomes()
        self.extract_genomes()
        self.genome_lineages()
        self.genome_metadata()
        self.build_taxdump()
        self.build_taxonmap()
        self.compile_database()

        if hasattr(self, 'mkdtemp'):
            rmtree(self.tmpdir)

        print(f'Database building finished at {timestamp()}.')

    def set_parameters(self, args):
        for key, val in vars(args).items():
            setattr(self, key, val)

        for key in 'capital', 'block', 'latin':
            get_config(self, key, f'taxonomy.{key}')
        for key in 'retries', 'delay', 'timeout':
            get_config(self, key, f'download.{key}')
        for key in 'diamond', 'makeblastdb':
            get_config(self, key, f'program.{key}')
        for key in 'threads', 'tmpdir':
            get_config(self, key, f'local.{key}')

        for key in 'capital', 'latin':
            setattr(self, key, arg2bool(getattr(self, key, None)))

        if not self.tmpdir:
            self.tmpdir = mkdtemp()
            setattr(self, 'mkdtemp', True)
        if not isdir(self.tmpdir):
            raise ValueError(f'Invalid temporary directory: {self.tmpdir}.')

        # Check executables (only for local commands, not for downloads)
        for key, exe in {'blast': 'makeblastdb', 'diamond': 'diamond'}.items():
            if self.compile in (key, 'both'):
                if getattr(self, exe) is None:
                    setattr(self, exe, exe)
                if which(getattr(self, exe)) is None:
                    raise ValueError(f'Invalid {exe} executable: {getattr(self, exe)}.')

        if self.default:
            print('The default protocol is selected for database building.')
            print('The program will download all protein sequences of NCBI RefSeq genomes of bacteria, archaea, fungi and protozoa, keep one genome per species, plus all NCBI-defined reference and representative genomes.')
            self.cats = 'microbe'
            self.sample = 1
            self.rank = 'species'
            self.reference = True
            self.represent = True

            if self.diamond is None:
                self.diamond = 'diamond'
            if which(self.diamond) is None:
                print('WARNING: Cannot find DIAMOND in this computer. You will need to manually compile the database after download is complete.')
                self.compile = 'none'
            else:
                self.compile = 'diamond'

        # Determine number of CPUs
        if self.compile in ('diamond', 'both') and not self.threads:
            try:
                self.threads = len(os.sched_getaffinity(0))
            except AttributeError:
                self.threads = os.cpu_count() or 1

        makedirs(self.output, exist_ok=True)
        self.down = join(self.output, 'download')
        makedirs(self.down, exist_ok=True)

    # ------------------------------------------------------------------------
    # Download helper using urllib (replaces rsync)
    # ------------------------------------------------------------------------
    def _download_file(self, url, dest_path, retries=None, delay=None):
        """Download a file from URL to dest_path with retries."""
        if retries is None:
            retries = self.retries if hasattr(self, 'retries') else 3
        if delay is None:
            delay = self.delay if hasattr(self, 'delay') else 5

        for attempt in range(retries):
            try:
                req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(req, timeout=self.timeout if hasattr(self, 'timeout') else 30) as resp:
                    with open(dest_path, 'wb') as f:
                        f.write(resp.read())
                return True
            except (urllib.error.URLError, urllib.error.HTTPError, OSError) as e:
                if attempt < retries - 1:
                    print(f"  Download failed: {e}. Retrying in {delay} seconds...")
                    sleep(delay)
                else:
                    print(f"  ERROR: Could not download {url} after {retries} attempts.")
                    raise
        return False

    # ------------------------------------------------------------------------
    # retrieve_taxdump: use HTTPS
    # ------------------------------------------------------------------------
    def retrieve_taxdump(self):
        fname = 'taxdump.tar.gz'
        url = f'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/{fname}'
        lfile = join(self.down, fname)

        if not self.check_local_file(lfile, self.overwrite):
            print('Downloading NCBI taxonomy database...', end='', flush=True)
            self._download_file(url, lfile)
            print(' done.')

        print('Reading NCBI taxonomy database...', end='', flush=True)
        with tarfile.open(lfile, 'r:gz') as f:
            f.extract('names.dmp', self.tmpdir)
            f.extract('nodes.dmp', self.tmpdir)
        self.taxdump = read_taxdump(self.tmpdir)
        print(' done.')
        print(f'  Total number of TaxIDs: {len(self.taxdump)}.')

    # ------------------------------------------------------------------------
    # retrieve_summary: use HTTPS
    # ------------------------------------------------------------------------
    def retrieve_summary(self, genbank=False):
        def get_summary(target):
            key = target.lower()
            fname = f'assembly_summary_{key}.txt'
            url = f'https://ftp.ncbi.nlm.nih.gov/genomes/{key}/{fname}'
            lfile = join(self.down, fname)

            if not self.check_local_file(lfile, self.overwrite):
                print(f'Downloading {target} assembly summary...', end='', flush=True)
                self._download_file(url, lfile)
                print(' done.')

            print(f'Reading {target} assembly summary...', end='', flush=True)
            df = pd.read_table(lfile, skiprows=1, low_memory=False)
            col = df.columns[0]
            if col.startswith('#'):
                df.rename(columns={col: col.lstrip('# ')}, inplace=True)
            print(' done.')
            return df

        self.df = get_summary('RefSeq')
        if self.genbank:
            self.df = pd.concat([self.df, get_summary('GenBank')])
        print(f'  Total number of genomes: {self.df.shape[0]}.')

    # ------------------------------------------------------------------------
    # retrieve_categories: remove rsync dependency; direct HTTPS download
    # ------------------------------------------------------------------------
    def retrieve_categories(self):
        cats = self.cats
        if cats == 'all':
            return
        elif cats == 'microbe':
            self.cats = ['archaea', 'bacteria', 'fungi', 'protozoa']
        else:
            self.cats = cats.split(',')
        print(f'Genome categories: {", ".join(self.cats)}')

        def get_categories(target):
            key = target.lower()
            makedirs(join(self.down, 'cats'), exist_ok=True)
            ldir = join(self.down, 'cats')
            fname = 'assembly_summary.txt'
            asms = []

            for cat in self.cats:
                lfile = join(ldir, f'{key}_{cat}.txt')
                url = f'https://ftp.ncbi.nlm.nih.gov/genomes/{key}/{cat}/{fname}'

                if self.check_local_file(lfile):
                    with open(lfile, 'r') as f:
                        asms_ = f.read().splitlines()
                else:
                    print(f'Downloading {target} {cat} assembly list...', flush=True)
                    try:
                        self._download_file(url, lfile)
                    except Exception:
                        raise ValueError(f'Category "{cat}" not found for {target} (invalid category or network issue).')
                    with open(lfile, 'r') as f:
                        asms_ = [x.split('\t', 1)[0] for x in f.read().splitlines() if not x.startswith('#')]

                print(f'  {cat}: {len(asms_)}')
                asms += asms_
            return asms

        asmset = set(get_categories('RefSeq'))
        if self.genbank:
            asmset.update(get_categories('GenBank'))
        self.df = self.df[self.df['assembly_accession'].isin(asmset)]
        print(f'  Total number of genomes in categories: {self.df.shape[0]}.')

    # ------------------------------------------------------------------------
    # filter_genomes, sort_genomes, filter_by_taxonomy, sample_by_taxonomy,
    # sample_by_quality, filter_to_sampled: unchanged from original (but keep)
    # ------------------------------------------------------------------------
    def filter_genomes(self):
        print('Filtering genomes...')
        n = self.df.shape[0]

        def report_diff(msg):
            nonlocal n
            n_ = self.df.shape[0]
            if n_ < n:
                print('  ' + msg.format(n - n_))
            n = n_

        if self.complete:
            self.df = self.df[self.df['assembly_level'].isin({'Complete Genome', 'Chromosome'})]
            report_diff('Dropped {} non-complete genomes.')

        self.df.rename(columns={'assembly_accession': 'accession'}, inplace=True)
        self.df['accnov'] = self.df['accession'].str.split('.', n=1).str[0]
        self.df['genome'] = 'G' + self.df['accnov'].str.split('_', n=1).str[-1]
        self.df.drop_duplicates(subset=['genome'], inplace=True)

        if self.genoids:
            self.genoids = set(list_from_param(self.genoids))
            print(f'{"Ex" if self.exclude else "In"}cluding {len(self.genoids)} custom genome IDs...')
            self.df = self.df[(self.df['accession'].isin(self.genoids) |
                               self.df['accnov'].isin(self.genoids) |
                               self.df['genome'].isin(self.genoids)) != self.exclude]
            report_diff('Dropped {} genomes.')

            if not self.exclude:
                self.capital = False
                self.latin = False
                self.block = ''
                self.sample = None

        self.df.query('ftp_path != "na"', inplace=True)
        report_diff('Dropped {} genomes without download link.')
        print('Done.')

    def sort_genomes(self):
        self.df['rc_seq'] = self.df.apply(
            lambda x: 0 if x['refseq_category'] == 'reference genome'
            else (1 if x['refseq_category'] == 'representative genome'
                  else (2 if pd.notnull(x['relation_to_type_material'])
                  else 3)), axis=1)

        self.df['al_seq'] = self.df['assembly_level'].map(
            {'Chromosome': 0, 'Complete Genome': 0, 'Scaffold': 1, 'Contig': 2})

        self.df.sort_values(by=['rc_seq', 'al_seq', 'genome'], inplace=True)
        self.df.drop(columns=['al_seq', 'rc_seq'], inplace=True)

    def filter_by_taxonomy(self):
        print('Filtering genomes by taxonomy...')
        n = self.df.shape[0]

        def report_diff(msg):
            nonlocal n
            n_ = self.df.shape[0]
            if n_ < n:
                print('  ' + msg.format(n - n_))
            n = n_

        if self.capital:
            self.df = self.df[self.df['organism_name'].apply(is_capital)]
            report_diff('Dropped {} genomes without capitalized organism name.')

        if self.block:
            self.block = list_from_param(self.block)
            self.df = self.df[~self.df['organism_name'].apply(contain_words, args=(self.block,))]
            report_diff('Dropped {} genomes with one or more blocked words in organism name.')

        self.df.drop(columns=['species_taxid'], inplace=True)

        self.df.dropna(subset=['taxid'], inplace=True)
        self.df['taxid'] = self.df['taxid'].astype(str)
        self.df = self.df[self.df['taxid'].isin(self.taxdump)]
        report_diff('Dropped {} genomes without valid taxId.')

        self.df['species'] = self.df['taxid'].apply(taxid_at_rank, rank='species', taxdump=self.taxdump)

        self.df.dropna(subset=['species'], inplace=True)
        report_diff('Dropped {} genomes without valid species taxId.')

        self.df['latin'] = self.df['species'].apply(lambda x: is_latin(self.taxdump[x]['name']))
        if self.latin:
            self.df.query('latin', inplace=True)
            report_diff('Dropped {} genomes without Latinate species name.')
        print('Done.')

        if self.taxids:
            self.taxids = set(list_from_param(self.taxids))
            print(f'{"Ex" if self.exclude else "In"}cluding {len(self.taxids)} custom TaxIDs...')
            self.df = self.df[self.df['taxid'].apply(lambda x: is_ancestral(x, self.taxids, self.taxdump) != self.exclude)]
            report_diff('Dropped {} genomes.')

    def sample_by_taxonomy(self):
        if self.sample is None:
            self.selected = set(self.df['genome'])
            return

        print('Sampling genomes based on taxonomy...')
        self.selected, n = set(), 0
        rank, sample, latin = self.rank, self.sample, False

        if rank == 'species_latin':
            rank, latin = 'species', True
        print(f'Up to {sample} genome(s) will be sampled per {rank}' + (' (Latinate names only) ' if latin else '') + '.')

        if rank not in self.df.columns:
            self.df[rank] = self.df['taxid'].apply(taxid_at_rank, rank=rank, taxdump=self.taxdump)
        df_ = self.df.dropna(subset=[rank])

        if latin:
            df_ = df_.query('latin')

        self.selected = set(df_.groupby(rank).head(sample)['genome'])
        n = df_[rank].unique().shape[0]
        print(f'  Sampled {len(self.selected)} genomes from {n} {rank_plural(rank)}.')

        if not self.above:
            return

        ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']
        if rank not in ranks:
            raise ValueError(f'Cannot determine taxonomic ranks above "{rank}", because it is not among the six standard ranks: {", ".join(ranks)}.')
        ranks = ranks[:ranks.index(rank)][::-1]
        print(f'Sampling will also be performed on: {", ".join(ranks)}.')

        for r in ranks:
            self.df['tmp_rank'] = self.df['taxid'].apply(taxid_at_rank, rank=r, taxdump=self.taxdump)
            df_ = self.df.dropna(subset=['tmp_rank'])

            goal = pd.Series(sample, df_['tmp_rank'].unique())
            done = df_.query('genome in @self.selected')['tmp_rank'].value_counts()
            left = goal.subtract(done, fill_value=0)
            left = left[left > 0].astype(int)

            s = []
            for tid, num in left.items():
                s.extend(df_.query(f'tmp_rank == "{tid}"').head(num)['genome'])
            print(f'  Sampled {len(s)} more genomes from {left.shape[0]} {rank_plural(r)}.')
            self.selected.update(s)

        print(f'  Sampled a total of {len(self.selected)} genomes at {rank} and above.')
        self.df.drop(columns=['tmp_rank'], inplace=True)

    def sample_by_quality(self):
        if self.reference:
            df_ = self.df.query('refseq_category == "reference genome"')
            print(f'Included {df_.shape[0]} reference genomes.')
            self.selected.update(df_['genome'])

        if self.represent:
            df_ = self.df.query('refseq_category == "representative genome"')
            print(f'Included {df_.shape[0]} representative genomes.')
            self.selected.update(df_['genome'])

        if self.typemater:
            df_ = self.df[self.df['relation_to_type_material'].notna()]
            print(f'Included {df_.shape[0]} type material genomes.')
            self.selected.update(df_['genome'])

    def filter_to_sampled(self):
        n = len(self.selected)
        if n == 0:
            raise ValueError('No genome is retained after sampling.')
        self.df.query('genome in @self.selected', inplace=True)
        self.df.sort_values('genome', inplace=True)
        print(f'Total number of sampled genomes: {n}.')

    # ------------------------------------------------------------------------
    # download_genomes: use HTTPS with urllib
    # ------------------------------------------------------------------------
    def download_genomes(self):
        if self.manual:
            fname = 'urls.txt'
            self.df['ftp_path'].to_csv(join(self.output, 'urls.txt'), header=False, index=False)
            print(f'URLs of sampled genomes written to {fname}. You may manually download their protein sequence data to faa/, then restart the database building pipeline.')
            sys.exit(0)

        print('Downloading non-redundant genomic data from NCBI...', flush=True)
        ldir = join(self.down, 'faa')
        makedirs(ldir, exist_ok=True)
        failed = []

        for row in self.df.itertuples():
            g = row.genome
            rdir = row.ftp_path.split('/', 3)[-1]
            stem = rdir.rsplit('/', 1)[-1]
            fname = f'{stem}_protein.faa.gz'
            lfile = join(ldir, fname)

            if self.check_local_file(lfile):
                continue

            # Build HTTPS URL from the ftp_path (replace ftp:// with https://)
            url = row.ftp_path.replace('ftp://', 'https://') + f'/{fname}'

            success = False
            for i in range(self.retries):
                try:
                    self._download_file(url, lfile, retries=1)  # single attempt inside loop
                    print('  ' + g, flush=True)
                    success = True
                    break
                except Exception:
                    if i < self.retries - 1:
                        sleep(self.delay)
                    else:
                        break

            if not success:
                print(f'WARNING: Cannot retrieve {fname}.')
                failed.append(g)

        print('Done.')

        if len(failed):
            print('Failed to retrieve the following genomes:')
            print('  ' + ', '.join(failed))
            failed_set = set(failed)
            self.df.query('genome not in @failed_set', inplace=True)

    # ------------------------------------------------------------------------
    # extract_genomes, genome_lineages, genome_metadata, build_taxdump,
    # build_taxonmap, compile_database: unchanged
    # ------------------------------------------------------------------------
    def extract_genomes(self):
        print('Extracting downloaded genomic data...', end='', flush=True)
        ldir = join(self.down, 'faa')
        prots = {}
        cp = None
        fout = open(join(self.output, 'db.faa'), 'w')

        def write_prot():
            try:
                fout.write(f'>{cp} {prots[cp]["name"]}\n{prots[cp]["seq"]}\n')
            except KeyError:
                return
            else:
                prots[cp]['aa'] = len(prots[cp]['seq'])
                del prots[cp]['name']
                del prots[cp]['seq']

        g2n, g2aa = {}, {}
        for row in self.df.itertuples():
            g, tid = row.genome, row.taxid
            g2n[g], g2aa[g] = 0, 0
            stem = row.ftp_path.rsplit('/', 1)[-1]
            lfile = join(ldir, f'{stem}_protein.faa.gz')
            with gzip.open(lfile, 'rb') as f:
                try:
                    content = f.read().decode().splitlines()
                except (OSError, EOFError, TypeError):
                    print(f' skipping corrupted file {stem}.', end='', flush=True)
                    continue
            cp = None
            for line in content:
                if line.startswith('>'):
                    write_prot()
                    p, name = line[1:].split(None, 1)
                    g2n[g] += 1
                    if p in prots:
                        cp = None
                        prots[p]['gs'].add(g)
                        prots[p]['tids'].add(tid)
                        g2aa[g] += prots[p]['aa']
                    else:
                        cp = p
                        prots[p] = {'name': name, 'gs': {g}, 'tids': {tid}, 'seq': ''}
                elif cp:
                    line = line.rstrip('*')
                    prots[cp]['seq'] += line
                    g2aa[g] += len(line)
            write_prot()
        fout.close()
        print(' done.')
        print('Combined protein sequences written to db.faa.')

        self.df['proteins'] = self.df['genome'].map(g2n)
        self.df['residues'] = self.df['genome'].map(g2aa)
        self.p2tids = {k: v['tids'] for k, v in prots.items()}

        print(f'  Total number of genomes extracted: {len(g2n)}.')
        print(f'  Total number of unique proteins extracted: {len(prots)}.')
        print('  Number of proteins shared by multiple genomes: {}.'.format(sum([1 for k, v in prots.items() if len(v['gs']) > 1])))
        print('  Number of proteins shared by multiple TaxIDs: {}.'.format(sum([1 for k, v in prots.items() if len(v['tids']) > 1])))

        fname = 'genome.map.gz'
        with gzip.open(join(self.output, fname), 'wb') as f:
            for k, v in sorted(prots.items()):
                f.write('{}\t{}\n'.format(k, ','.join(sorted(v['gs']))).encode())
        print(f'Protein-to-genome(s) map written to {fname}.')

    def genome_lineages(self):
        ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        self.df[ranks[:-1]] = self.df['taxid'].apply(lambda x: pd.Series(taxids_at_ranks(x, ranks[:-1], self.taxdump)))

        print('Number of taxonomic groups represented:')
        for rank in ranks:
            print(f'  {rank}: {self.df[rank].nunique()}.')

        self.df['kingdom'] = self.df[['superkingdom', 'kingdom']].apply(lambda x: x[1] if x[1] else x[0], axis=1)

        self.df['lineage'] = self.df[ranks[1:]].fillna('').apply(lambda col: col.apply(lambda val, name: '{}__{}'.format(name[0], self.taxdump[val]['name'] if val else ''), args=(col.name,))).apply('; '.join, axis=1)

        fname = 'lineages.txt'
        self.df['lineage'].to_csv(join(self.output, fname), sep='\t', header=False)
        print(f'Genome lineages written to {fname}.')

    def genome_metadata(self):
        self.df.set_index('genome', inplace=True)
        self.df = self.df[['proteins', 'residues', 'assembly_level', 'accession', 'bioproject', 'biosample', 'asm_name', 'organism_name', 'infraspecific_name', 'isolate', 'taxid', 'ftp_path']]
        fname = 'genomes.tsv'
        self.df.to_csv(join(self.output, fname), sep='\t')
        print(f'Genome metadata written to {fname}.')

    def build_taxdump(self):
        tids = self.df['taxid'].unique().tolist()

        ancs = set(tids)
        for tid in tids:
            cid = tid
            while True:
                pid = self.taxdump[cid]['parent']
                if cid == pid or pid in ancs:
                    break
                ancs.add(pid)
                cid = pid

        dir_ = join(self.output, 'taxdump')
        makedirs(dir_, exist_ok=True)
        for fname in ('nodes.dmp', 'names.dmp'):
            fo = open(join(dir_, fname), 'w')
            fi = open(join(self.tmpdir, fname), 'r')
            for line in fi:
                row = line.rstrip('\r\n').replace('\t|', '').split('\t')
                if row[0] in ancs:
                    if fname == 'nodes.dmp' or 'scientific name' in row:
                        fo.write(line)
            fi.close()
            fo.close()
        print('Taxonomy database written to taxdump/.')

    def build_taxonmap(self):
        self.taxonmap = {p: max(tids) if len(tids) == 1 else find_lca(tids, self.taxdump) for p, tids in self.p2tids.items()}

        fname = 'taxon.map.gz'
        with gzip.open(join(self.output, fname), 'wb') as f:
            for p, tid in sorted(self.taxonmap.items()):
                f.write(f'{p}\t{tid}\n'.encode())
        print(f'Protein-to-taxonomy map written to {fname}.')

    def compile_database(self):
        if self.compile == 'none':
            return
        if self.compile in ('blast', 'both'):
            self.build_blast_db()
        if self.compile in ('diamond', 'both'):
            self.build_diamond_db()

    def build_blast_db(self):
        taxonmap = join(self.tmpdir, 'taxon.map')
        with open(taxonmap, 'w') as f:
            for p, tid in sorted(self.taxonmap.items()):
                f.write(f'{p}\t{tid}\n')

        makedirs(join(self.output, 'blast'), exist_ok=True)
        cmd = ' '.join((
            self.makeblastdb,
            '-dbtype', 'prot',
            '-in', join(self.output, 'db.faa'),
            '-out', join(self.output, 'blast', 'db'),
            '-title', 'db',
            '-parse_seqids',
            '-taxid_map', taxonmap))
        print('Build BLAST database...', flush=True)
        ec, out = run_command(cmd)
        if ec:
            raise ValueError(f'makeblastdb failed with error code {ec}.')

        remove(taxonmap)
        print('Done.')

    def build_diamond_db(self):
        taxonmap = join(self.tmpdir, 'prot.accession2taxid')
        with open(taxonmap, 'w') as f:
            f.write('accession\taccession.version\ttaxid\tgi\n')
            for p, tid in sorted(self.taxonmap.items()):
                f.write(f'{p.rsplit(".", 1)[0]}\t{p}\t{tid}\tna\n')

        makedirs(join(self.output, 'diamond'), exist_ok=True)
        cmd = ' '.join((
            self.diamond, 'makedb',
            '--threads', str(self.threads),
            '--in', join(self.output, 'db.faa'),
            '--taxonmap', taxonmap,
            '--taxonnodes', join(self.output, 'taxdump', 'nodes.dmp'),
            '--taxonnames', join(self.output, 'taxdump', 'names.dmp'),
            '--db', join(self.output, 'diamond', 'db')))
        print('Build DIAMOND database...', flush=True)
        ec, out = run_command(cmd)
        if ec:
            raise ValueError(f'diamond failed with error code {ec}.')

        remove(taxonmap)
        print('Done.')

    @staticmethod
    def check_local_file(file, overwrite=False):
        if isfile(file):
            if getsize(file) == 0:
                remove(file)
            elif overwrite:
                print(f'  WARNING: existing local file {basename(file)} will be overwritten.')
                remove(file)
            else:
                print(f'  Using local file {basename(file)}.')
                return True
        return False

# ----------------------------------------------------------------------------
# Command-line entry point (if run standalone)
# ----------------------------------------------------------------------------
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    # Build argument parser from the 'arguments' structure
    for group in arguments:
        if isinstance(group, str):
            continue
        for arg in group:
            if isinstance(arg, list):
                flags = arg[0].split('|')
                kwargs = arg[2]
                parser.add_argument(*flags, help=arg[1], **kwargs)
    args = parser.parse_args()
    db = Database()
    db(args)
```

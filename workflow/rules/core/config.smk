from collections import OrderedDict

import contextlib
import typing
import types
from snakemake.utils import logger, validate
from snakemake.io import _load_configfile
from dataclasses import dataclass, field, asdict

WORKFLOW_DIR = str(workflow.current_basedir)
SCHEMA_DIR = os.path.realpath(
    os.path.join(WORKFLOW_DIR, os.pardir, os.pardir, "schemas")
)


class PropertyDict(OrderedDict):
    """Simple class that allows for property access"""

    def __init__(self, data=dict()):
        super().__init__(data)
        if isinstance(data, types.GeneratorType):
            return
        for k, v in data.items():
            if isinstance(v, dict):
                v = PropertyDict(v)
            elif isinstance(v, list):
                val = []
                for x in v:
                    if isinstance(x, PropertyDict):
                        val.append(x)
                    elif isinstance(x, dict):
                        val.append(PropertyDict(x))
                    else:
                        val.append(x)
                v = val
            else:
                pass
            self[k] = v

    def __setitem__(self, key, value):
        OrderedDict.__setitem__(self, key, value)
        if key not in dir(dict()):
            try:
                setattr(self, key, value)
            except Exception as e:
                print(e)
                print(key, value)
                raise


class Schema(PropertyDict):
    def __init__(self, schemafile):
        self.schemafile = schemafile
        data = _load_configfile(self.schemafile, filetype="Schema")
        super().__init__(data)

    def validate(self, data):
        validate(data, schema=self.schemafile)


assemblies_schema = Schema(os.path.join(SCHEMA_DIR, "assemblies.schema.yaml"))
reads_schema = Schema(os.path.join(SCHEMA_DIR, "reads.schema.yaml"))
transcripts_schema = Schema(os.path.join(SCHEMA_DIR, "transcripts.schema.yaml"))
genecovr_schema = Schema(os.path.join(SCHEMA_DIR, "genecovr_csv.schema.yaml"))
definitions_schema = Schema(os.path.join(SCHEMA_DIR, "definitions.schema.yaml"))
tools_schema = Schema(os.path.join(SCHEMA_DIR, "tools.schema.yaml"))
config_schema = Schema(os.path.join(SCHEMA_DIR, "config.schema.yaml"))


class SampleData:
    _index = ["id"]
    _idcols = ["id"]
    schema = None

    def __init__(self, *args):
        if len(args) == 0:
            self._data = pd.DataFrame(columns=self._index)
        elif len(args) == 1:
            args = args[0]
            if isinstance(args, SampleData):
                self._index = args._index
                self.schema = args.schema
                self._data = args.data
            else:
                self._load_data(args)
        elif len(args) > 1:
            assert all(isinstance(x, SampleData) for x in args), logger.error(
                "all instances must be SampleData"
            )
            self._index = args[0]._index
            self.schema = args[0].schema
            self._data = pd.concat(x.data for x in args)
        else:
            raise TypeError
        self.schema.validate(self.data)
        if "id" not in self._data.columns and self._index == ["id"]:
            logger.info(f"generating id column from {self._idcols}")
            self._data["id"] = "_".join(self._data[self._idcols])
            self._data["id"] = self._data[self._idcols].agg("_".join, axis=1)
        self._data.set_index(self._index, inplace=True, drop=False)
        self._data = self._data.replace({np.nan: None})
        self._data.index.names = self._index

    def _load_data(self, infile):
        if infile is None:
            self._data = pd.DataFrame(columns=self._index)
            return
        ext = os.path.splitext(infile)[1]
        if ext == ".yaml":
            with open(infile) as fh:
                data = yaml.load(fh, yaml.Loader)
            assert isinstance(data, list)
            self._data = pd.DataFrame(data)
        elif ext == ".tsv":
            self._data = pd.read_csv(infile, sep="\t")

    def subset(self, invert=False, **kw):
        keep = self.data.index.isin(self.data.index)
        for k, v in kw.items():
            if v is None:
                continue
            if isinstance(v, set):
                v = list(v)
            if not isinstance(v, list):
                v = [v]
            if len(v) == 0:
                continue
            try:
                if k in self.data.index.names:
                    keep = keep & self.data.index.isin(v)
                else:
                    keep = keep & self.data[k].isin(v)
            except KeyError as e:
                print(e)
                raise
        cls = type(self)
        new = cls(self)
        if invert:
            keep = ~keep
        new._data = new._data[keep]
        return new

    @property
    def data(self):
        return self._data

    @property
    def ids(self):
        return self.data.index.to_list()


class Assemblies(SampleData):
    _idcols = ["species", "version"]
    schema = assemblies_schema


class Reads(SampleData):
    schema = reads_schema


class Transcripts(SampleData):
    schema = transcripts_schema


@dataclass
class Tool:
    _default = None
    _all_input = None
    _analysis = None
    name: str

    def __post_init__(self):
        self._fmt = {
            "input": self._all_input,
        }

    def _extra(self):
        return

    def to_dict(self):
        d = {**asdict(self), **self._analysis.to_dict()}
        if self._extra() is not None:
            d.update(**self._extra())
        return d

    def format(self, key):
        fmt = self._fmt[key]
        return expand(fmt, **self.to_dict())

    @property
    def input(self) -> list[str]:
        return self.format("input")


@dataclass
class Blast(Tool):
    _default = tools_schema.definitions["blast.config"].properties
    database: list[str]


@dataclass
class Busco(Tool):
    _default = tools_schema.definitions["busco.config"].properties
    _all_input = "{results}/busco/{analysis}/{assembly_ids}/{mode}/run_{lineage}/short_summary_{assembly_ids}.txt"
    lineage: str
    mode: str = _default.mode.default


@dataclass
class Jellyfish(Tool):
    _default = tools_schema.definitions["jellyfish.config"].properties
    _all_input = "{results}/jellyfish/{analysis}/{assembly_ids}/merged.{kmer}_jf.hist"
    tmpdir: str = None
    kmer: list[int] = field(default_factory=lambda: _default.kmer.default)
    npartitions: int = _default.npartitions.default
    count_pairs: str = _default.count_pairs.default

    def _extra(self):
        d = {"partition": range(0, self.npartitions)}
        return d

    def format(self, key):
        fmt = self._fmt[key]
        val = expand(fmt, **self.to_dict())
        if self.count_pairs:
            pfx = "{results}/jellyfish/{analysis}/kmer_comparison/{assembly_ids}.{kmer}_jf.png"
            val.extend(expand(pfx, **self.to_dict()))
        return val

    def count_input(self, wildcards):
        if wildcards.dataset in self._analysis.assembly_ids:
            return {
                "seq": "{interim}/jellyfish/{analysis}/{dataset}/{prefix}.fasta".format(
                    **wildcards
                )
            }
        return {"seq": self._analysis.get_reads()}

    def merge_input(self, wildcards):
        if wildcards.dataset in self._analysis.assembly_ids:
            fmt = "{{interim}}/jellyfish/{analysis}/{dataset}/{{partition}}.{kmer}mer_counts.jf".format(
                **wildcards
            )
            return {"jf": expand(fmt, **self.to_dict())}
        fmt = "{{interim}}/jellyfish/{analysis}/{dataset}/{{prefix}}.{kmer}mer_counts.jf".format(
            **wildcards
        )
        return {
            "jf": expand(
                fmt,
                prefix=[os.path.basename(x) for x in self._analysis.get_reads()],
                **self.to_dict(),
            )
        }


@dataclass
class Kraken2(Tool):
    _default = tools_schema.definitions["kraken2.config"].properties
    _all_input = (
        "{results}/kraken2/{analysis}/{assembly_ids}/{db}.{window_size}.report.txt"
    )

    db: str
    window_size: list[int] = field(default_factory=lambda: _default.window_size.default)
    npartitions: int = _default.npartitions.default

    def _extra(self):
        d = {"partition": range(0, self.npartitions)}
        return d

    def results(self, wildcards):
        fmt = f"{{interim}}/kraken2/{wildcards.analysis}/{wildcards.assembly}/{wildcards.db}.{wildcards.length}.{{partition}}.{{suffix}}"
        d = {
            "output": expand(fmt, suffix="output.txt.gz", **self.to_dict()),
            "unclassified": expand(
                fmt, suffix="unclassified.fasta.gz", **self.to_dict()
            ),
        }
        return d

    def reports(self, wildcards):
        fmt = f"{{interim}}/kraken2/{wildcards.analysis}/{wildcards.assembly}/{wildcards.db}.{wildcards.length}.{{partition}}.report.txt"
        return expand(fmt, **self.to_dict())


@dataclass
class Quast(Tool):
    _default = None
    _all_input = "{results}/quast/{analysis}/{assembly_ids}/{rpt}"

    rpt: list[str] = field(
        default_factory=lambda: [
            "report.tsv",
            "transposed_report.tsv",
            "report.txt",
            "transposed_report.txt",
        ]
    )


@dataclass
class Genecovr(Tool):
    _default = tools_schema.definitions["genecovr.config"].patternProperties.properties
    _all_input = Path("{results}/genecovr/{analysis}/")
    __GENECOVR_MINMATCH__ = [0.75, 0.85, 0.9, 0.95]
    __GENECOVR_NCONTIGS__ = ["bar"]
    __GENECOVR_MATCH_INDEL__ = ["violin", "boxplot", "boxplot.log10"]
    __GENECOVR_FN__ = ["width_violin.pdf", "qnuminsert.pdf"]
    __GENECOVR_DEPTH_BREADTH__ = ["coverage", "jitter", "hist", "seqlengths"]
    __GENECOVR_CSV_GZ__ = [
        "gene_body_coverage.csv.gz",
        "psldata.csv.gz",
        "gbc_summary.csv.gz",
        "ncontigs_per_transcripts.csv.gz",
    ]

    csvfile: str = _default.csvfile.default
    outprefix: str = _default.outprefix.default

    def __post_init__(self):
        super().__post_init__()
        if self.csvfile is None:
            self.csvfile = self._all_input / "genecovr.csv"

    def format(self, key):
        pfx = self._fmt[key]
        d = self.to_dict()
        retval = []
        retval += report(
            expand(
                pfx / "gene_body_coverage.minmatch.{mm}.pdf",
                mm=self.__GENECOVR_MINMATCH__,
                **d,
            ),
            caption="../report/genecovr_gbc.rst",
            category="Gene body coverages",
        )
        retval += report(
            expand(
                pfx / "ncontigs_per_transcripts.{type}.mm0.75.pdf",
                type=self.__GENECOVR_NCONTIGS__,
                **d,
            ),
            caption="../report/genecovr_ncontigs.rst",
            category="Number of contigs per transcript",
        )
        retval += report(
            expand(
                pfx / "depth_breadth_{type}.mm0.75.pdf",
                type=self.__GENECOVR_DEPTH_BREADTH__,
                **d,
            ),
            caption="../report/genecovr_depth_breadth.rst",
            category="Depth and breadth of coverage",
        )
        retval += report(
            expand(
                pfx / "match_indel.{type}.pdf", type=self.__GENECOVR_MATCH_INDEL__, **d
            ),
            caption="../report/genecovr_match_indel.rst",
            category="Match and indel distributions",
        )
        retval += report(
            expand(pfx / "{fn}", fn=self.__GENECOVR_FN__, **d),
            caption="../report/genecovr_match_indel.rst",
            category="Match and indel distributions",
        )
        retval += report(
            expand(pfx / "{fn}", fn=self.__GENECOVR_CSV_GZ__, **d),
            caption="../report/genecovr_data.rst",
            category="Data files",
        )
        return retval

    def csv_input(self):
        interim = self.to_dict()["interim"]
        assembly_ids = self._analysis.assemblies.ids
        transcript_ids = self._analysis.transcripts.ids
        assembly_fasta = self._analysis.assemblies.data["fasta"].tolist()
        trx_fasta = self._analysis.transcripts.data["fasta"].tolist()
        retval = {
            "dataset": [
                f"{a}/{b}" for a, b in itertools.product(assembly_ids, transcript_ids)
            ],
            "psl": [
                str(f"{interim}/gmap/map/{a}-{b}.psl")
                for a, b in itertools.product(assembly_ids, transcript_ids)
            ],
            "assembly": [
                f"{a}.fai" for a, b in itertools.product(assembly_fasta, trx_fasta)
            ],
            "trxset": [f"{b}" for a, b in itertools.product(assembly_fasta, trx_fasta)],
        }
        return retval


@dataclass
class Repeatmasker(Tool):
    _default = tools_schema.definitions["repeatmasker.config"].properties


tools = {
    "blast": Blast,
    "busco": Busco,
    "jellyfish": Jellyfish,
    "kraken2": Kraken2,
    "quast": Quast,
    "genecovr": Genecovr,
}


@dataclass
class Analysis:
    name: str = None
    label: str = None
    description: str = None
    fs: dict = field(default_factory=PropertyDict)
    assembly_ids: list = None
    read_ids: list = None
    transcript_ids: list = None
    assemblies: list = field(default_factory=Assemblies)
    reads: list = field(default_factory=Reads)
    transcripts: list = field(default_factory=Transcripts)
    tools: dict = field(default_factory=PropertyDict)

    _section = "analysis"

    def __post_init__(self):
        # Need to update assembly_ids etc if empty list, no?
        self.assemblies = self.assemblies.subset(id=self.assembly_ids)
        self.reads = self.reads.subset(id=self.read_ids)
        self.transcripts = self.transcripts.subset(id=self.transcript_ids)
        self.name = self.name.replace(self._section + os.sep, "")
        for k, v in self.tools.items():
            if isinstance(v, bool):
                v = {}
            t = tools[k](name=k, **v)
            t._analysis = self
            self.tools[k] = t

    @property
    def analysis_tools(self):
        return self.tools.keys()

    def to_dict(self):
        d = {
            "analysis": self.name,
            "assembly_ids": self.assemblies.ids,
            "read_ids": self.reads.ids,
            "transcript_ids": self.transcripts.ids,
            **self.fs,
        }
        return d

    def get_reads(self):
        """Retrieve the sequence files for a set of read ids"""
        return (
            self.reads.data["read1"].to_list()
            + self.reads.data["read2"].dropna().to_list()
        )

    def get_assembly(self, assembly_id, fai=False):
        fasta = self.assemblies.subset(id=assembly_id).data.fasta.to_list()
        assert len(fasta) == 1
        if fai:
            return "{}.fai".format(fasta[0])
        return fasta[0]

    def get_transcriptome(self, transcript_id):
        fasta = self.transcripts.subset(id=assembly_id).data.fasta.to_list()
        assert len(fasta) == 1
        return fasta[0]


@dataclass
class RuleConfig:
    _default = workflow.default_resources.parsed
    name: str
    envmodules: list[str] = field(default_factory=list)
    options: typing.Union[str, list, dict] = ""
    runtime: int = _default.get("runtime", 100)
    threads: int = 1
    attempt: int = 1
    mem_mb: int = _default.get("mem_mb", 1000)
    java_options: str = _default.get("java_options", "")
    java_tmpdir: str = _default.get("java_tmpdir", _default["tmpdir"])
    window_size: list[int] = field(default_factory=list)
    step_size: list[int] = field(default_factory=list)

    def xthreads(self, wildcards, attempt):
        return attempt * self.threads

    def xruntime(self, wildcards, attempt):
        return attempt * self.runtime

    def xmem(self, wildcards, attempt):
        return attempt * self.mem_mb


class Config(PropertyDict):
    _analysissection = "analysis"
    _assemblies: Assemblies
    _reads: Reads
    _transcripts: Transcripts

    def __init__(self, conf, *args, **kw):
        super().__init__(conf)
        self.__post_init__()

    def __post_init__(self):
        self._assemblies = Assemblies(config["assemblies"])
        self._transcripts = Transcripts(config["transcripts"])
        self._reads = Reads(config["reads"])
        self._init_analyses()

    def _init_analyses(self):
        regex = self._analysissection + os.sep
        for k in self.keys():
            if not k.startswith(regex):
                continue
            d = {"name": k, **self[k], **self.to_dict()}
            self[k] = Analysis(**d)
        # Add default analysis with tools
        self["analysis/default"] = Analysis(
            **{"name": "analysis/default", **self.default, **self.to_dict()}
        )

    def ruleconf(self, rulename, **kw):
        """Retrieve rule configuration"""
        data = {"name": rulename, **kw}
        if "rules" in self.keys():
            if rulename in self.rules:
                data.update(**self.rules[rulename])
        return RuleConfig(**data)

    @property
    def analyses(self):
        return [v for k, v in self.items() if isinstance(v, Analysis)]

    def analyses_w_tool(self, tool):
        return [a for a in self.analyses if tool in a.analysis_tools]

    @property
    def analysisnames(self):
        return [x.name for x in self.analyses]

    def analysis(self, name):
        return self[f"{self._analysissection}/{name}"]

    @property
    def default(self):
        d = {}
        akeys = Analysis.__dataclass_fields__.keys()
        for k in self.keys():
            if k in akeys:
                d.update(**{k: self[k]})
        return d

    def to_dict(self):
        return {
            "reads": self._reads,
            "assemblies": self._assemblies,
            "transcripts": self._transcripts,
            "fs": self.fs,
        }

    def get_assembly(self, assembly_id, fai=False):
        fasta = self._assemblies.subset(id=assembly_id).data.fasta.to_list()
        assert len(fasta) == 1
        if fai:
            return "{}.fai".format(fasta[0])
        return fasta[0]

    def get_transcriptome(self, transcript_id):
        fasta = self._transcripts.subset(id=transcript_id).data.fasta.to_list()
        assert len(fasta) == 1
        return fasta


# context manager for cd
@contextlib.contextmanager
def cd(path, logger):
    CWD = os.getcwd()
    logger.info("Changing directory from {} to {}".format(CWD, path))

    os.chdir(path)
    try:
        yield
    except Exception as e:
        logger.warning(e)
        logger.warning("Exception caught: ".format(sys.exc_info()[0]))
    finally:
        logger.info("Changing directory back to {}".format(CWD))
        os.chdir(CWD)

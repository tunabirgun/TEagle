import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "app" / "web" / "app.js"


def _source():
    return APP.read_text(encoding="utf-8")


def _export_helpers():
    source = _source()
    start = source.index("function sha256Hex")
    end = source.index("window.exportTable", start)
    return source[start:end]


def _node(expression):
    prelude = """
global.window={}; global.state={}; global.seqEl={value:""}; global.err=()=>{};
global.downloadBlob=()=>{}; global.Blob=class Blob{};
"""
    script = prelude + _export_helpers() + "\n" + expression
    result = subprocess.run(["node", "-e", script], cwd=ROOT, text=True,
                            capture_output=True, check=False)
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_web_export_helpers_are_valid_javascript():
    result = subprocess.run(["node", "--check", str(APP)], cwd=ROOT, text=True,
                            capture_output=True, check=False)
    assert result.returncode == 0, result.stderr


def test_model_gff3_bed_fasta_coordinate_contract_and_order():
    observed = _node(r"""
const model={exons:[{start:12,end:20,strand:"+"},{start:0,end:6,strand:"+"}],
 introns:[{start:6,end:12,strand:"+",donor:"GT",acceptor:"AG",canonical:true}],
 cds:[{start:2,end:6,strand:"+",phase:0,cds_group:"cds_1"},{start:12,end:18,strand:"+",phase:2,cds_group:"cds_1"}]};
const seq="ACGT".repeat(5), g=modelGff3(model,seq,"chr1","fetched_gene_model"),
 b=modelBed(model,seq.length,"chr1"), f=modelFasta({exons:[{start:4,end:8,strand:"-"}]},"AAAACCGT","chr1",["exon"]);
console.log(JSON.stringify({
 g:g.split("\n").filter(x=>x&&!x.startsWith("#")&&x.includes("\t")).map(x=>{const p=x.split("\t");return p.slice(2,7).concat(p[7]);}),
 b:b.trim().split("\n").map(x=>x.split("\t").slice(1,3)), fasta:f.trim().split("\n")[1]
}));
""")
    assert observed["g"] == [["exon", "1", "6", ".", "+", "."],
                              ["CDS", "3", "6", ".", "+", "0"],
                              ["intron", "7", "12", ".", "+", "."],
                              ["CDS", "13", "18", ".", "+", "2"],
                              ["exon", "13", "20", ".", "+", "."]]
    assert observed["b"] == [["0", "6"], ["2", "6"], ["6", "12"],
                              ["12", "18"], ["12", "20"]]
    assert observed["fasta"] == "ACGG"


def test_model_export_mutations_fail_and_path_like_id_is_not_a_path():
    observed = _node(r"""
const failures=[];
for(const f of [{start:-1,end:2,strand:"+"},{start:0,end:0,strand:"+"},
 {start:0,end:9,strand:"+"},{start:0,end:2,strand:"?"}]){
 try{modelBed({exons:[f]},8,"../../sample A;unit=1\\x");failures.push(false);}catch(e){failures.push(true);}
}
const id="../../sample A;unit=1\\x", model={exons:[{start:0,end:2,strand:"+"}]};
const g=modelGff3(model,"ACGT",id,"fetched_gene_model"), b=modelBed(model,4,id);
console.log(JSON.stringify({failures,id:safeRecordId(id),stem:safeFileStem(id),g,b}));
""")
    assert observed["failures"] == [True, True, True, True]
    assert "/" not in observed["id"] and "\\" not in observed["id"]
    assert "/" not in observed["stem"] and "\\" not in observed["stem"]
    import hashlib
    digest = hashlib.sha256("../../sample A;unit=1\\x".encode()).hexdigest()
    assert f"original_record_sha256={digest}" in observed["g"]
    assert f"record_sha256={digest}" in observed["b"]
    assert "../../sample" not in observed["g"] and "../../sample" not in observed["b"]


def test_web_cds_gff3_fails_closed_without_phase_but_bed_and_fasta_remain():
    observed = _node(r"""
const model={cds:[{start:0,end:3,strand:"+"}]};
let gffFailed=false; try{modelGff3(model,"ACGT","chr1","fetched_gene_model");}catch(e){gffFailed=/phase/.test(e.message);}
    console.log(JSON.stringify({gffFailed,bed:modelBed(model,4,"chr1"),fasta:modelFasta(model,"ACGT","chr1",["CDS"])}));
""")
    assert observed["gffFailed"] is True
    import hashlib
    digest = hashlib.sha256(b"chr1").hexdigest()
    assert "\tCDS_1|record_sha256=" in observed["bed"]
    assert f">record_{digest[:12]}.cds1 record_sha256={digest}" in observed["fasta"]


def test_web_never_reemits_windows_unix_or_traversal_record_ids():
    observed = _node(r"""
const ids=["C:\\Users\\person\\private\\sample.fa","C:private_sample.fa","/home/person/private/sample.fa","../../private/sample.fa"];
const model={exons:[{start:0,end:2,strand:"+"}]};
console.log(JSON.stringify(ids.map(id=>({id,safe:safeRecordId(id),stem:safeFileStem(id),
 g:modelGff3(model,"ACGT",id,"fetched_gene_model"),b:modelBed(model,4,id),f:modelFasta(model,"ACGT",id,["exon"]),digest:sha256Hex(id)}))));
""")
    for item in observed:
        assert item["id"] not in item["g"] + item["b"] + item["f"]
        assert item["digest"] in item["g"] and item["digest"] in item["b"] and item["digest"] in item["f"]
        assert item["safe"] == "record_" + item["digest"][:12]
        assert item["stem"] == "record_" + item["digest"]


def test_web_local_identity_policy_hashes_every_header_even_when_it_looks_public():
    observed = _node(r"""
const publicIds=["M11240","M11240.1","AM040263","NC_000013.11","GCF_009914755.1",
 "chr1","chr1:100-200","input_sequence","record_deadbeefdead","P1"];
const sensitive=["patient_042","donor-A1","sample_tuna_birgun","descriptive specimen header","patient_042 (input)"];
const model={exons:[{start:0,end:2,strand:"+"}]};
console.log(JSON.stringify({
 local:publicIds.concat(sensitive).map(raw=>({raw,identity:safeExportIdentity(raw),stem:safeFileStem(raw),
  g:modelGff3(model,"ACGT",raw,"fetched_gene_model"),b:modelBed(model,4,raw),f:modelFasta(model,"ACGT",raw,["exon"])}))
}));
""")
    for item in observed["local"]:
        import hashlib
        digest = hashlib.sha256(item["raw"].encode()).hexdigest()
        assert item["identity"] == {
            "id": f"record_{digest[:12]}", "idSha256": digest,
            "internalKey": f"sha256:{digest}", "label": f"record_{digest[:12]}", "role": "",
        }
        assert item["stem"] == f"record_{digest}"
        combined = item["g"] + item["b"] + item["f"]
        assert item["raw"] not in combined
        assert digest in item["g"] and digest in item["b"] and digest in item["f"]


def test_web_accepts_real_backend_compact_full_digest_shape_without_rehashing():
    import hashlib
    digest = hashlib.sha256(b"patient_042").hexdigest()
    compact = f"record_{digest[:12]}"
    observed = _node(f"""
const identity=safeExportIdentity({json.dumps(compact)},{json.dumps(digest)});
const role=safeExportIdentity({json.dumps(compact + ' (input)')},{json.dumps(digest)});
const model={{exons:[{{start:0,end:2,strand:'+'}}]}};
console.log(JSON.stringify({{identity,role,stem:safeFileStem({json.dumps(compact)},'record',{json.dumps(digest)}),
 g:modelGff3(model,'ACGT',{json.dumps(compact)},'fetched_gene_model',{json.dumps(digest)}),
 b:modelBed(model,4,{json.dumps(compact)},{json.dumps(digest)}),
 f:modelFasta(model,'ACGT',{json.dumps(compact)},['exon'],{json.dumps(digest)})}}));
""")
    assert observed["identity"] == {
        "id": compact, "label": compact, "role": "", "idSha256": digest,
        "internalKey": f"sha256:{digest}",
    }
    assert observed["role"]["label"] == compact + " (input)"
    assert observed["role"]["role"] == "input"
    assert observed["stem"] == f"record_{digest}"
    assert digest in observed["g"] and digest in observed["b"] and digest in observed["f"]


def test_web_full_digest_binding_disambiguates_a_compact_id_collision():
    observed = _node(r"""
const a="aaaaaaaaaaaa"+"1".repeat(52),b="aaaaaaaaaaaa"+"2".repeat(52),compact="record_aaaaaaaaaaaa";
console.log(JSON.stringify({identities:[safeExportIdentity(compact,a),safeExportIdentity(compact+" (custom background)",b)],
 stems:[safeFileStem(compact,"record",a),safeFileStem(compact+" (custom background)","record",b)]}));
""")
    identities = observed["identities"]
    assert identities[0]["id"] == identities[1]["id"] == "record_aaaaaaaaaaaa"
    assert identities[0]["idSha256"] != identities[1]["idSha256"]
    assert identities[0]["internalKey"] != identities[1]["internalKey"]
    assert observed["stems"] == ["record_" + "a" * 12 + "1" * 52,
                                  "record_" + "a" * 12 + "2" * 52]


def test_web_rejects_mismatched_or_malformed_compact_digest_join():
    observed = _node(r"""
const failures=[];
for(const [id,digest] of [["record_aaaaaaaaaaaa","b".repeat(64)],["record_aaaaaaaaaaa","a".repeat(64)],
 ["NC_000013.11","a".repeat(64)],["record_aaaaaaaaaaaa","A".repeat(64)]]){
 try{safeExportIdentity(id,digest);failures.push(false);}catch(e){failures.push(/does not match/.test(e.message));}
}
console.log(JSON.stringify(failures));
""")
    assert observed == [True, True, True, True]


def test_generated_identity_requires_bare_input_and_a_validated_result_binding():
    observed = _node(r"""
const attempts=[];
for(const submitted of [">input_sequence\nACGT",">P1\nACGT",">record_deadbeefdead\nACGT"]){
 const id=submitted.slice(1).split("\n",1)[0],record={id,identitySource:"generated",seq:"ACGT",composition:{length:4}};
 try{validateAnalysisRecordIdentities({records:[record]},submitted,null);attempts.push(false);}catch(e){attempts.push(/generated/.test(e.message));}
}
const loose=safeExportIdentity("input_sequence");
console.log(JSON.stringify({attempts,loose}));
""")
    assert observed["attempts"] == [True, True, True]
    assert observed["loose"]["id"].startswith("record_")
    assert observed["loose"]["id"] != "input_sequence"


def test_web_canonical_source_url_accepts_only_accession_bound_ncbi_or_ena_url():
    observed = _node(r"""
const cases=[
 ["https://www.ncbi.nlm.nih.gov/nuccore/NC_000013.11","NC_000013.11","NCBI nuccore"],
 ["https://www.ebi.ac.uk/ena/browser/view/AM040263.1","AM040263.1","ENA fallback"],
 ["https://www.ncbi.nlm.nih.gov/nuccore/NC_000014.9","NC_000013.11","NCBI nuccore"],
 ["https://www.ncbi.nlm.nih.gov/nuccore/NC_000013.11?report=fasta","NC_000013.11","NCBI nuccore"],
 ["https://www.ncbi.nlm.nih.gov/nuccore/NC_000013.11#x","NC_000013.11","NCBI nuccore"],
 ["https://www.ncbi.nlm.nih.gov/nuccore/NC_000013.11.evil","NC_000013.11","NCBI nuccore"],
 ["https://www.ncbi.nlm.nih.gov/nuccore/NC_000013.11/extra","NC_000013.11","NCBI nuccore"],
 ["https://www.ebi.ac.uk/ena/browser/view/NC_000013.11","NC_000013.11","NCBI nuccore"],
 ["https://evil.example/NC_000013.11","NC_000013.11","NCBI nuccore"],
 ["javascript:alert(1)","NC_000013.11","NCBI nuccore"],
 ["https://www.ncbi.nlm.nih.gov/nuccore/patient_042","patient_042","NCBI nuccore"]
];
console.log(JSON.stringify(cases.map(x=>canonicalSourceUrl(...x))));
""")
    assert observed[:2] == [
        "https://www.ncbi.nlm.nih.gov/nuccore/NC_000013.11",
        "https://www.ebi.ac.uk/ena/browser/view/AM040263.1",
    ]
    assert observed[2:] == [""] * 9


def test_fetched_accession_trust_requires_complete_schema_and_cannot_be_forged():
    observed = _node(r"""
const payload={accession:"NC_000013.11",source:"NCBI nuccore",
 endpoint:"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
 sourceUrl:"https://www.ncbi.nlm.nih.gov/nuccore/NC_000013.11",
 retrievedUtc:"2026-08-22T01:02:03+00:00",taxid:9606,fasta:">NC_000013.11 chromosome 13\nACGT\n"};
payload.sourceContext={sourceContextVersion:"1.0",inputType:"fetched_accession",source:"ncbi_nuccore",endpoint:"ncbi_efetch",accession:payload.accession,retrievedUtc:payload.retrievedUtc,sourceUrl:payload.sourceUrl,records:[{id:payload.accession,length:4,sha256:sha256Hex("ACGT")}],taxid:"9606"};
const context=validatedFetchedSource(payload), forged={...context};
const mutations=[
 {...payload,accession:"NC_000014.9"},{...payload,source:"NCBI"},
 {...payload,endpoint:payload.endpoint+"?id=NC_000013.11"},
 {...payload,sourceUrl:payload.sourceUrl+"?report=fasta"},
 {...payload,sourceUrl:payload.sourceUrl+"#fragment"},
 {...payload,sourceUrl:payload.sourceUrl+".evil"},
 {...payload,retrievedUtc:"2026-08-22"},{...payload,taxid:"9606x"},
 {...payload,taxid:0},{...payload,length:5},
 {...payload,fasta:">NC_000013.11\nAC\n>NC_000013.11\nGT\n"},
 {...payload,fasta:">patient_042\nACGT\n"}
];
console.log(JSON.stringify({context,frozen:Object.isFrozen(context),keys:Object.keys(context).sort(),trusted:safeExportIdentity(payload.accession,null,context,"ACGT"),
 trustedStem:safeFileStem(payload.accession,"record",null,context,"ACGT"),
 local:safeExportIdentity(payload.accession),forged:safeExportIdentity(payload.accession,null,forged,"ACGT"),
 wrongSequence:safeExportIdentity(payload.accession,null,context,"TGCA"),
 rejected:mutations.map(validatedFetchedSource)}));
""")
    assert observed["context"] == {
        "sourceContextVersion": "1.0", "inputType": "fetched_accession", "source": "ncbi_nuccore",
        "endpoint": "ncbi_efetch", "accession": "NC_000013.11", "taxid": "9606",
        "retrievedUtc": "2026-08-22T01:02:03+00:00",
        "sourceUrl": "https://www.ncbi.nlm.nih.gov/nuccore/NC_000013.11",
        "records": [{
            "id": "NC_000013.11", "length": 4,
            "sha256": "1dff3e84fe7877e0673b69bbddcf40124e396e3f9943dd890c91b6a09adb9af0",
        }],
    }
    assert observed["keys"] == sorted(observed["context"])
    assert observed["frozen"] is True
    assert observed["trusted"]["id"] == "NC_000013.11"
    assert "idSha256" not in observed["trusted"]
    assert observed["trustedStem"] == "NC_000013.11"
    assert observed["local"]["id"].startswith("record_")
    assert observed["forged"]["id"].startswith("record_")
    assert observed["wrongSequence"]["id"].startswith("record_")
    assert observed["rejected"] == [None] * 12


def test_ena_fetch_context_requires_exact_endpoint_and_emits_enum_only():
    observed = _node(r"""
const payload={accession:"AM040263.1",source:"ENA (EMBL-EBI)",
 endpoint:"https://www.ebi.ac.uk/ena/browser/api/fasta/AM040263.1",
 sourceUrl:"https://www.ebi.ac.uk/ena/browser/view/AM040263.1",
 retrievedUtc:"2026-08-22T01:02:03Z",taxid:"3702",fasta:">AM040263.1 description\nACGT\n"};
payload.sourceContext={sourceContextVersion:"1.0",inputType:"fetched_accession",source:"ena_fasta",endpoint:"ena_fasta",accession:payload.accession,retrievedUtc:payload.retrievedUtc,sourceUrl:payload.sourceUrl,records:[{id:payload.accession,length:4,sha256:sha256Hex("ACGT")}],taxid:"3702"};
const context=validatedFetchedSource(payload);
console.log(JSON.stringify({context,identity:safeExportIdentity(payload.accession,null,context,"ACGT"),
 bad:validatedFetchedSource({...payload,endpoint:"https://www.ebi.ac.uk/ena/browser/api/fasta/AM040263.1?download=true"})}));
""")
    assert observed["context"]["source"] == "ena_fasta"
    assert observed["context"]["endpoint"] == "ena_fasta"
    assert observed["identity"]["id"] == "AM040263.1"
    assert observed["bad"] is None


def test_source_wire_is_exact_serializable_v1_and_forged_variants_downgrade_to_local():
    observed = _node(r"""
const raw={accession:"NC_000013.11",source:"NCBI nuccore",
 endpoint:"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
 sourceUrl:"https://www.ncbi.nlm.nih.gov/nuccore/NC_000013.11",
 retrievedUtc:"2026-08-22T01:02:03Z",taxid:9606,fasta:">NC_000013.11 record\nACGT\n"};
raw.sourceContext={sourceContextVersion:"1.0",inputType:"fetched_accession",source:"ncbi_nuccore",endpoint:"ncbi_efetch",accession:raw.accession,retrievedUtc:raw.retrievedUtc,sourceUrl:raw.sourceUrl,records:[{id:raw.accession,length:4,sha256:sha256Hex("ACGT")}],taxid:"9606"};
const context=validatedFetchedSource(raw),wire=sourceWire(context);
const variants=[
 {...wire},{...wire,records:[]},{...wire,records:[{...wire.records[0],length:5}]},
 {...wire,records:[{...wire.records[0],sha256:"0".repeat(64)}]},
 {...wire,accession:"NC_000014.9"},{...wire,sourceContextVersion:"0.9"},
 {...wire,extra:"forged"},{inputType:"fetched_accession",accession:raw.accession}
];
console.log(JSON.stringify({wire,wireJson:JSON.stringify(wire),local:sourceWire(null),
 downgraded:variants.map(sourceWire),rawStored:wire!==raw}));
""")
    wire = observed["wire"]
    assert set(wire) == {
        "sourceContextVersion", "inputType", "source", "endpoint", "accession",
        "retrievedUtc", "sourceUrl", "records", "taxid",
    }
    assert wire["records"] == [{
        "id": "NC_000013.11", "length": 4,
        "sha256": "1dff3e84fe7877e0673b69bbddcf40124e396e3f9943dd890c91b6a09adb9af0",
    }]
    assert json.loads(observed["wireJson"]) == wire
    local = {"sourceContextVersion": "1.0", "inputType": "local_input"}
    assert observed["local"] == local
    assert observed["downgraded"] == [local] * 8
    assert observed["rawStored"] is True
    from teagle_core import provenance
    records = [("NC_000013.11 record", "ACGT")]
    assert provenance.validate_source_context(wire, records) == wire
    assert provenance.validate_source_context({**wire, "extra": "forged"}, records) == local


def test_coordinate_source_wire_matches_backend_schema_and_rejects_forged_bindings():
    observed = _node(r"""
const source={accession:"NC_000013.11",organism:"patient_042",taxid:"9606",
 assemblyAccession:"GCF_000001405.40",assemblyName:"donor-A1",
 regions:[{chrAccession:"NC_000013.11",start:100,stop:103,strand:1},
          {chrAccession:"NC_000014.9",start:200,stop:203,strand:2}],
 coordSystem:"1-based-inclusive",retrievalType:"coordinate",displayLocus:"private locus",
 source:"NCBI E-utilities (efetch seq_start/seq_stop)",
 endpoint:"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
 retrievedUtc:"2026-08-22T10:20:30Z"};
const raw={ok:true,fasta:">NC_000013.11:100-103\nACGT\n>NC_000014.9:200-203(-)\nTGCA\n",
 seq_length:8,organism:"patient_042",assemblyName:"donor-A1",source};
raw.sourceContext={sourceContextVersion:"1.0",inputType:"fetched_coordinate",source:"ncbi_coordinate",
 endpoint:"ncbi_efetch_range",accession:"NC_000013.11",assemblyAccession:"GCF_000001405.40",
 retrievedUtc:"2026-08-22T10:20:30Z",coordSystem:"1-based-inclusive",
 regions:[{chrAccession:"NC_000013.11",start:100,end:103,strand:1},
          {chrAccession:"NC_000014.9",start:200,end:203,strand:2}],
 sourceUrl:"https://www.ncbi.nlm.nih.gov/nuccore/NC_000013.11",taxid:"9606",
 records:[{id:"NC_000013.11:100-103",length:4,sha256:sha256Hex("ACGT")},
          {id:"NC_000014.9:200-203",length:4,sha256:sha256Hex("TGCA")} ]};
const context=validatedFetchedSource(raw),wire=sourceWire(context),mutations=[
 {...raw.sourceContext,extra:"x"},
 {...raw.sourceContext,sourceUrl:"https://www.ncbi.nlm.nih.gov/nuccore/NC_000013.11?report=fasta"},
 {...raw.sourceContext,records:[raw.sourceContext.records[1],raw.sourceContext.records[0]]},
 {...raw.sourceContext,records:[raw.sourceContext.records[0],{...raw.sourceContext.records[1],sha256:"0".repeat(64)}]},
 {...raw.sourceContext,regions:[raw.sourceContext.regions[1],raw.sourceContext.regions[0]]},
 {...raw.sourceContext,taxid:9606}
];
console.log(JSON.stringify({wire,identities:[
 safeExportIdentity("NC_000013.11:100-103",null,context,"ACGT"),
 safeExportIdentity("NC_000014.9:200-203",null,context,"TGCA")],
 rejected:mutations.map(sourceContext=>validatedFetchedSource({...raw,sourceContext})),
 rawFields:JSON.stringify(wire)}));
""")
    wire = observed["wire"]
    assert set(wire) == {
        "sourceContextVersion", "inputType", "source", "endpoint", "accession",
        "assemblyAccession", "retrievedUtc", "coordSystem", "regions", "sourceUrl",
        "records", "taxid",
    }
    assert wire["regions"] == [
        {"chrAccession": "NC_000013.11", "start": 100, "end": 103, "strand": 1},
        {"chrAccession": "NC_000014.9", "start": 200, "end": 203, "strand": 2},
    ]
    assert [identity["id"] for identity in observed["identities"]] == [
        "NC_000013.11:100-103", "NC_000014.9:200-203",
    ]
    assert observed["rejected"] == [None] * 6
    for private_value in ("patient_042", "donor-A1", "private locus", "(-)"):
        assert private_value not in observed["rawFields"]

    from teagle_core import provenance
    records = [("NC_000013.11:100-103", "ACGT"), ("NC_000014.9:200-203(-)", "TGCA")]
    assert provenance.validate_source_context(wire, records) == wire
    local = {"sourceContextVersion": "1.0", "inputType": "local_input"}
    backend_mutations = [
        {**wire, "extra": "x"},
        {**wire, "sourceUrl": wire["sourceUrl"] + "?report=fasta"},
        {**wire, "taxid": 9606},
        {**wire, "records": list(reversed(wire["records"]))},
    ]
    assert [provenance.validate_source_context(value, records) for value in backend_mutations] == [local] * 4


def test_every_web_operation_sends_the_normalized_source_wire():
    source = _source()
    assert 'source: state.source' not in source and 'source:state.source' not in source
    for route in ("/api/analyze", "/api/annotate", "/api/splice", "/api/primers", "/api/pcr"):
        assert route in source
    assert source.count("source: sourceWire()") >= 5
    assert "r.useSource ? sourceWire() : LOCAL_SOURCE_WIRE" in source
    assert "const rawFetch = await api(\"/api/fetch\"" in source


def test_web_exposes_every_analysis_record_and_fails_closed_for_unverified_json():
    source = _source()
    assert "d.records.map((r,i)=>" in source
    assert "All records are listed" in source
    assert "complete result JSON unavailable" in source
    assert "cannot independently verify the backend seal" in source
    assert 'downloadText(strictJson(state.lastAnalysis)' not in source
    assert 'downloadText(strictJson(state.lastProvenance)' not in source
    assert "state.recordIndex!=null&&state.recordIndex!==index" in source
    assert "clearRecordBoundOutputs" in source
    assert "records[0]" not in source
    assert "first is used downstream" not in source
    assert "first record annotated" not in source


def test_stale_and_unverified_export_clicks_fail_without_downloading():
    helpers = _export_helpers()
    script = """
global.window={};global.state={};global.seqEl={value:""};global.errors=[];global.downloads=[];
global.err=x=>errors.push(x);global.downloadBlob=(blob,name)=>downloads.push({parts:blob.parts,name});
global.Blob=class Blob{constructor(parts){this.parts=parts;}};
""" + helpers + r"""
window.exportSelectedRecordFasta(); window.exportGeneModel("gff3");
state.lastAnalysis={resultSchemaVersion:"1.0"};state.lastProvenance={manifestSha256:"0".repeat(64)};
window.exportAnalysisJson();window.exportProvenanceJson();
console.log(JSON.stringify({errors,downloads}));
"""
    result = subprocess.run(["node", "-e", script], cwd=ROOT, text=True,
                            capture_output=True, check=False)
    assert result.returncode == 0, result.stderr
    observed = json.loads(result.stdout)
    assert len(observed["errors"]) == 4
    assert observed["downloads"] == []


def test_selected_record_fasta_requires_engine_bound_identity_and_supports_all_trust_classes():
    helpers = _export_helpers()
    script = """
global.window={};global.state={};global.seqEl={value:""};
global.errors=[];global.downloads=[];global.err=x=>errors.push(x);
global.downloadBlob=(blob,name)=>downloads.push({text:blob.parts.join(""),name});
global.Blob=class Blob{constructor(parts){this.parts=parts;}};
""" + helpers + r"""
for(const id of ["C:\\Users\\person\\private.fa","patient_042","donor-A1","sample_tuna_birgun","descriptive specimen header"]){
 state.lastRec={id,seq:"ACGT"};window.exportSelectedRecordFasta();
}
const raw="actual backend local header",digest=sha256Hex(raw),id=`record_${digest.slice(0,12)}`;
let rec={id,idSha256:digest,identitySource:"local_input",seq:"ACGT",composition:{length:4}};
validateAnalysisRecordIdentities({records:[rec]},`>${raw}\nACGT`,null);state.lastRec=rec;window.exportSelectedRecordFasta();
rec={id:"input_sequence",identitySource:"generated",seq:"ACGT",composition:{length:4}};
validateAnalysisRecordIdentities({records:[rec]},"ACGT",null);state.lastRec=rec;window.exportSelectedRecordFasta();
const fetched={accession:"NC_000013.11",source:"NCBI nuccore",endpoint:"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",sourceUrl:"https://www.ncbi.nlm.nih.gov/nuccore/NC_000013.11",retrievedUtc:"2026-08-22T01:02:03Z",taxid:9606,fasta:">NC_000013.11 description\nACGT\n"};
fetched.sourceContext={sourceContextVersion:"1.0",inputType:"fetched_accession",source:"ncbi_nuccore",endpoint:"ncbi_efetch",accession:fetched.accession,retrievedUtc:fetched.retrievedUtc,sourceUrl:fetched.sourceUrl,records:[{id:fetched.accession,length:4,sha256:sha256Hex("ACGT")}],taxid:"9606"};
const context=validatedFetchedSource(fetched);
rec={id:"NC_000013.11",identitySource:"validated_fetch",seq:"ACGT",composition:{length:4}};
validateAnalysisRecordIdentities({records:[rec]},fetched.fasta,context);state.source=context;state.lastRec=rec;window.exportSelectedRecordFasta();
console.log(JSON.stringify({errors,downloads}));
"""
    result = subprocess.run(["node", "-e", script], cwd=ROOT, text=True,
                            capture_output=True, check=False)
    assert result.returncode == 0, result.stderr
    observed = json.loads(result.stdout)
    assert len(observed["errors"]) == 5 and len(observed["downloads"]) == 3
    import hashlib
    digest = hashlib.sha256(b"actual backend local header").hexdigest()
    backend_export, generated_export, fetched_export = observed["downloads"]
    assert backend_export["name"] == f"record_{digest}.fasta"
    assert f">record_{digest[:12]} original_record_sha256={digest}" in backend_export["text"]
    assert generated_export == {"text": ">input_sequence\nACGT\n", "name": "input_sequence.fasta"}
    assert fetched_export == {"text": ">NC_000013.11\nACGT\n", "name": "NC_000013.11.fasta"}


def test_pcr_table_and_fasta_exports_bind_sensitive_sources_by_full_digest():
    source = _source()
    table = source[source.index("window.exportTable"):source.index("function csvBtns")]
    amplicons = source[source.index("window.copyAmplicons"):source.index("/* ---------- health")]
    script = ("global.window={};global.state={};global.seqEl={value:''};global._LC=new Set();"
              "global.familyCategory=()=>'';global.ampCall=()=>'<b>intended</b>';"
              "global.captured={};global.downloadTable=(headers,rows)=>captured.table={headers,rows};"
              "global.copyText=text=>captured.fasta=text;global.event=null;\n" +
              _export_helpers() + table + amplicons +
              "\nconst raws=['patient_042','donor-A1','sample_tuna_birgun','descriptive specimen header'];"
              "state={lastPcr:{amplicons:raws.map((raw,i)=>{const sourceIdSha256=sha256Hex(raw),id='record_'+sourceIdSha256.slice(0,12);return {source:id+(i?' (custom background)':' (input)'),sourceIdSha256,pair:'P'+(i+1),start:1,end:5,length:4,seq:'ACGT'};})}};"
              "window.exportTable('pcr','csv');window.copyAmplicons();console.log(JSON.stringify(captured));")
    result = subprocess.run(["node", "-e", script], cwd=ROOT, text=True,
                            capture_output=True, check=False)
    assert result.returncode == 0, result.stderr
    observed = json.loads(result.stdout)
    import hashlib
    raw_ids = ["patient_042", "donor-A1", "sample_tuna_birgun",
               "descriptive specimen header"]
    for raw, row in zip(raw_ids, observed["table"]["rows"]):
        digest = hashlib.sha256(raw.encode()).hexdigest()
        role = " (input)" if raw == raw_ids[0] else " (custom background)"
        assert row[3:5] == [f"record_{digest[:12]}{role}", digest]
        assert raw not in observed["fasta"]
        role_field = "input" if raw == raw_ids[0] else "custom_background"
        assert f"amplicon_record_{digest[:12]}_1-5_4bp source_role={role_field} source_sha256={digest}" in observed["fasta"]


def test_stale_csv_tsv_clicks_have_no_download_callback():
    source = _source()
    block = source[source.index("window.exportTable"):source.index("function csvBtns")]
    script = ("global.window={};global.state={lastRec:null,orfs:[],domains:[],family:[],candidates:[],lastPcr:null,features:null,splice:null};"
              "global.downloads=[];global.downloadTable=()=>downloads.push(1);global._LC=new Set();"
              "global.familyCategory=()=>\"other_or_unclassified\";\n" + block +
              "\nfor(const k of ['struct','orf','domain','family','primers','pcr','gene','splice'])window.exportTable(k,'csv');"
              "console.log(JSON.stringify(downloads));")
    result = subprocess.run(["node", "-e", script], cwd=ROOT, text=True,
                            capture_output=True, check=False)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == []


def test_all_mutation_paths_clear_state_dom_and_export_callbacks():
    source = _source()
    for text in ("Sequence changed — run analysis", "Fetch started — prior results",
                 "Upload started — prior results", "A new analysis is running",
                 "Analysis failed — no result", "Record selection changed"):
        assert text in source
    assert 'for(const id of ["structBody","wslBody","spliceBody","primBody","pcrBody","provBody"])' in source
    assert 'for(const id of ["figStruct","figGene","figSplice","figGel"])' in source
    for field in ("state.structural=[]", "state.orfs=[]", "state.domains=[]", "state.family=[]",
                  "state.lastRec=null", "state.candidates=[]", "state.lastPcr=null"):
        assert field in source


def test_invalidation_removes_live_export_dom_and_figure_callbacks():
    source = _source()
    block = source[source.index("function clearDisplayedResults"):source.index("function countNt")]
    script = r"""
const ids=["structBody","wslBody","spliceBody","primBody","pcrBody","provBody","figStruct","figGene","figSplice","figGel","structMeta","wslMeta","rRecords","rStruct","rOrf","mLen","mGC","mN","mValid","design","runpcr"];
const elements=Object.fromEntries(ids.map(id=>[id,{innerHTML:'<button onclick="download()">old export</button>',textContent:'old',disabled:false}]));
global.$=s=>elements[s.slice(1)]||null;global.esc=x=>x;global.renderPcrQueue=()=>{};
global._gvState={figStruct:{download:1},figGene:{download:1},figSplice:{download:1},figGel:{download:1}};
global._figState={figStruct:1,figGene:1,figSplice:1,figGel:1};
global.state={structural:[1],orfs:[1],domains:[1],family:[1],lastRec:{},candidates:[1],lastPcr:{},lastAnalysis:{},lastProvenance:{}};
""" + block + r"""
invalidateAllResults("stale");
console.log(JSON.stringify({state,html:Object.fromEntries(ids.slice(0,6).map(id=>[id,elements[id].innerHTML])),gv:_gvState,fig:_figState}));
"""
    result = subprocess.run(["node", "-e", script], cwd=ROOT, text=True,
                            capture_output=True, check=False)
    assert result.returncode == 0, result.stderr
    observed = json.loads(result.stdout)
    assert all("old export" not in html and "onclick" not in html for html in observed["html"].values())
    assert observed["gv"] == {} and observed["fig"] == {}
    assert observed["state"]["structural"] == [] and observed["state"]["lastRec"] is None
    assert observed["state"]["lastAnalysis"] is None and observed["state"]["lastProvenance"] is None


def test_family_rows_are_not_silently_dropped_or_mislabeled():
    source = _source()
    assert "state.family = hits" in source
    assert "const source=state.familySequence" in source
    assert "const rows = hits.map" in source
    assert "const f = state.family || []" in source
    assert "Every returned RepeatMasker row is retained" in source
    assert "an unclassified row is not evidence that the interval is non-TE" in source
    assert "Dfam 4.0 curated" not in source
    assert "Family annotation accepts exactly one FASTA record per run" in source
    assert '"classification_category"' in source and '"scope_note"' in source
    assert source.count('"scope_note"') >= 7
    assert 'function familyCategory(classFamily)' in source
    assert 'const rows = hits.map' in source and 'rows = f.map' in source
    assert 'familyCategory(x.class_family)' in source


def test_fetched_source_url_and_all_pair_manifests_are_retained():
    source = _source()
    assert "validatedFetchedSource(rawFetch)" in source
    assert "state.source=sourceContext" in source
    assert "sourceUrl:d.sourceUrl" not in source
    assert "renderProvenanceSet(manifests)" in source
    assert "Every staged pair returned its own manifest; all seals are listed below." in source
    assert "The manifest set is incomplete and export is refused." in source
    assert "showing P1" not in source


def test_truthful_web_export_surface_has_only_implemented_formats():
    source = _source()
    for token in ("exportGeneModel('fasta','exon')", "exportGeneModel('fasta','CDS')",
                  "exportGeneModel('gff3')", "exportGeneModel('bed')",
                  "exportSpliceModel('fasta','exon')", "exportSpliceModel('gff3')",
                  "exportSpliceModel('bed')"):
        assert token in source
    export_section = source[source.index("/* ---------- deterministic scientific text exports ---------- */"):
                            source.index("window.exportTable")]
    assert "xlsx" not in export_section.lower() and ".pdf" not in export_section.lower()

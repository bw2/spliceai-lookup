import json
import os
import re
from flask import Flask, request, Response
from flask_cors import CORS

from spliceai.utils import Annotator, get_delta_scores

ANNOTATOR = {
    "grch37": Annotator(os.path.expanduser("~/p1/ref/GRCh37/hg19.fa"), "grch37"),
    "grch38": Annotator(os.path.expanduser("~/p1/ref/GRCh38/hg38.fa"), "grch38"),
}

DISTANCE = 50  # maximum distance between the variant and gained/lost splice site, defaults to 50
MASK = 0  # mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss, defaults to 0

app = Flask(__name__, template_folder='.')
CORS(app)

VARIANT_RE = re.compile(
    "(chr)?(?P<chrom>[0-9XYMTt]{1,2})"
    "[: -]+"
    "(?P<pos>[0-9]{1,9})"
    "[: -]+"
    "(?P<ref>[ACGT]+)"
    "[: ->]+"
    "(?P<alt>[ACGT]+)"
)


def parse_variant(variant_str):
    match = VARIANT_RE.match(variant_str)
    if not match:
        raise ValueError(f"Unable to parse variant: {variant_str}")

    return match['chrom'], int(match['pos']), match['ref'], match['alt']


class VariantRecord:
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alts = [alt]

    def __repr__(self):
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"


API_ROUTE = "/splice-ai/<genome_version>"
EXAMPLE = f"For example: {API_ROUTE.replace('<genome_version>', 'grch37')}?variants='chr1:12345 A>G',chr2-12346-C-TG"

@app.route(API_ROUTE, methods=['POST', 'GET'])
def get_spliceai_scores(genome_version):
    if genome_version not in ANNOTATOR.keys():
        return f"Invalid genome version: {genome_version}. It should be {' or '.join(ANNOTATOR.keys())}\n", 400

    # check params
    params = {}
    if request.values:
        params.update(request.values)

    if 'variants' not in params:
        params.update(request.get_json(force=True, silent=True) or {})

    variants = params.get('variants')
    if not variants:
        return f'"variants" arg not specified. This api expects a "variants" arg as either a GET or a POST parameter. {EXAMPLE}\n', 400

    if isinstance(variants, str):
        variants = variants.split(",")

    if not isinstance(variants, list):
        return f'"variants" arg must be a comma-separated string or a json list. Instead found {type(variants)}: {variants}\n', 400

    # parse and perform liftover
    results = []
    for variant in variants:
        variant = variant.strip().strip("'").strip('"').strip(",")
        if not variant:
            continue

        try:
            chrom, pos, ref, alt = parse_variant(variant)
        except ValueError as e:
            results.append(f"ERROR: {e}")
            continue
            #return f'Unable to parse "{locus}": {str(e)}\n', 400

        record = VariantRecord(f"chr{chrom}", pos, ref, alt)
        scores = get_delta_scores(record, ANNOTATOR[genome_version], DISTANCE, MASK)
        print(scores)
        if len(scores) == 0:
            results.append(f"ERROR: unable to compute scores for {variant}")
        else:
            results.append(scores)

    return Response(json.dumps(results),  mimetype='application/json')


@app.route('/', defaults={'path': ''})
@app.route('/<path:path>')
def catch_all(path):

    return f"""<html>
<head>
<title>SpliceAI Lookup API</title>
</head>
<body style="font-family: monospace">
Welcome to the SpliceAI Lookup API. <br/>
<br />
The API can be queried via: <br />
<br />
GET  {API_ROUTE}?variants=[variant1],[variant2]... <br />
<br/>
or <br/>
<br/>
POST  /{API_ROUTE}<br />
{{variants: "[variant1],[variant2],[variant3]"}} <br/>
<br/>
{EXAMPLE} <br />
<br />
<br />
<b>[genome_version]</b> should be: {' or '.join(ANNOTATOR.keys())} <br />
<b>variants should have the format "chrom:pos ref&gt;alt" or "chrom-pos-ref-alt" or "chrom pos ref alt" <br />
<br />

The API response is a json list that's the same length as the input list and has splice AI scores for each variant.<br/>
<br />
</body>
</html>"""


if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))

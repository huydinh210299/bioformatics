from flask import Flask, request
from flask_restplus import Resource, Api, reqparse, fields, marshal
from flask_cors import CORS
from ML.classify import Classify
from werkzeug.exceptions import BadRequest
from werkzeug.datastructures import FileStorage
import json
import os


app = Flask(__name__)
CORS(app)
api = Api(app)

classification = Classify()

@api.route('/hello')
class HelloWorld(Resource):
    def get(Self):
        return {'hello': 'world'}

### /result?filename=xxx
@api.route('/result')
class Classify(Resource):
    def get(Self):
        args = request.args
        filename = args["filename"]
        check_fn = os.path.exists(filename)
        while(not check_fn):
            check_fn = os.path.exists(filename)
        try:
            rs =  classification.antiMicroable('./'+ filename)
            result = rs.to_json(orient="records")
            parsed = json.loads(result)
            json.dumps(parsed, indent=4)
            if os.path.exists(filename):
                os.remove(filename)
            filename = filename.split(".")[0] + '.fasta'
            if os.path.exists(filename):
                os.remove(filename)
            return parsed
        except:
            raise BadRequest()

### Thuc hien upload input la sequence
parser = reqparse.RequestParser()
parser.add_argument('seq', type=str)

@api.route('/upload-data/seq')
@api.expect(parser)
class UploadSeq(Resource):
    def post(Self):
        try:
            args = parser.parse_args()
            proseq = args['seq']
            fn_fasta = classification.paacBySeq(proseq)
            fn_csv = fn_fasta.strip('.fasta') + '.paac.csv'
            return {'filename': fn_csv}
        except:
            raise BadRequest()

### Thuc hien upload la file fasta
upload_parser = api.parser()
upload_parser.add_argument('file', location='files',
                           type=FileStorage, required=True)

@api.route('/upload-data/fastafile')
@api.expect(upload_parser)
class Upload(Resource):
    def post(self):
        try:
            args = upload_parser.parse_args()
            uploaded_file = args['file']
            fn = uploaded_file.filename
            uploaded_file.save(fn)
            fn_fasta = classification.paacByFile(fn)
            fn_csv = fn_fasta.split(".")[0] + '.paac.csv'
            return {'filename': fn_csv}
        except:
            raise BadRequest()
if __name__ == '__main__':
    app.run(debug=True)
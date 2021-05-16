
import { Button, Carousel, Input, Radio, Form, Upload, notification } from 'antd';
import axios from 'axios';
import './App.css';
import image1 from './image/image1.png';
import image2 from './image/image2.jpg';
import image3 from './image/image3.png';
import image4 from './image/image4.png';
import {CheckCircleOutlined, InboxOutlined, LoadingOutlined} from '@ant-design/icons';
import { useState } from 'react';

const validateMessages = {
  required: '${label} is required!',
  types: {
    email: '${label} is not a valid email!',
    number: '${label} is not a valid number!',
  },
  number: {
    range: '${label} must be between ${min} and ${max}',
  },
};


function App() {
  const [radioValue, setRadioValue] = useState('sequence');
  const [inputFile, setInputFile] = useState();
  const [requireFile, setRequireFile] = useState(false);
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState([]);


  const onChangeRadio = (e) => {
    console.log('radio', e.target.value);
    setRadioValue(e.target.value);
  };


  const onChangeFile = (e) => {
    console.log('files ',e.target.files);
    const files = e.target.files;
    const filesArr = Array.prototype.slice.call(files);
    setInputFile(filesArr[0]);
    setRequireFile(false);
  }

  const onFinish = async (values) => {
    setLoading(true);
    let res1;

    if(radioValue === 'sequence') {
      console.log('values', values.inputText);
      res1 = await axios.post('http://127.0.0.1:5000/upload-data/seq', {
        seq: values.inputText
      });
      console.log('resp1', res1.data.filename);
      
    } else {
      if(!inputFile){
        setLoading(false);
        setRequireFile(true);
        return;
      }else{
        setRequireFile(false);
        const formData = new FormData();
        formData.append('file', inputFile);
        res1 = await axios.post('http://127.0.0.1:5000/upload-data/fastafile', formData);
        console.log('resp1', res1.data.filename);
      }
    }

    setTimeout(async ()=>{
      setLoading(false);
      axios.get(`http://127.0.0.1:5000/result?filename=${res1.data.filename}`)
      .then((res) => {
        setResult(res.data);
        console.log('resp2', res.data);
      })
      .catch((err)=>{
        notification['error']({
          message: 'Request to server error!',
          description: 'Time out of server'
        });
      });
      
    }, 4000);
  }
  return (
    <div className="App">
      <Carousel autoplay>
        <div className="carousel-item">
          <img className="image" src={image2} alt="" />
        </div>
        <div className="carousel-item">
          <img className="image" src={image2} alt="" />
        </div>
        <div className="carousel-item">
          <img className="image" src={image3} alt="" />
        </div>
        <div className="carousel-item">
          <img className="image" src={image4} alt="" />
        </div>
      </Carousel>
      <h1>A web server to classify the antimicrobial protein and peptide</h1>
      <div>
        <Form name="nest-messages" onFinish={onFinish} validateMessages={validateMessages}>
          <h2>Input your sequense:</h2>
          {/* <div className="form-input">
            <textarea className="text-input" name="input" cols="60" rows="5" placeholder="Example: ACSAG"></textarea>
          </div> */}
          {
            radioValue === 'sequence' && (<Form.Item 
              name={['inputText']}
              rules={[
                {
                  required: true,
                },
              ]}
              label="Input text"
            >
              <Input.TextArea 
                placeholder="Example: ACSAG"
                style={{width: '400px', height: '150px', borderRadius: '10px', fontSize: '18px'}}
              />
            </Form.Item>)
          }
          {
            radioValue === 'fasta' && (
              <div
              style={{display: 'flex',flexDirection: 'column' , justifyContent: 'center', alignItems: 'center '}}
              >
                {/* <Form.Item >
                  <Form.Item 
                    valuePropName="fileList" 
                    name="dragger"
                    noStyle
                    label="Input file"
                    rules={[
                      {
                        required: true,
                      },
                    ]}
                  > */}
                    <label 
                      htmlFor="validate_other_dragger"
                      style={{width: '300px'}}
                    >
                      <input 
                        type="file" 
                        name="uploadFile"
                        onChange={onChangeFile}
                        id="validate_other_dragger"
                        style={{display: 'none'}}
                      />
                      <div
                        style={{background: '#fff', padding: '10px 20px', borderRadius: '10px', cursor: 'pointer'}}
                      >
                        <p
                          style={{color: 'blue', fontSize: '30px', marginBottom: '0'}}
                        >
                          <InboxOutlined />
                        </p>
                        <p style={{fontSize: '18px'}}>Click to this area chose file fasta to upload!</p>
                      </div>
                    </label>
                    {inputFile && (
                      <div style={{color:'#fff', fontSize: '20px'}}>{inputFile.name}</div>
                    )}
                    {requireFile && (
                      <div style={{color:'red', fontSize: '20px'}}>Input file is require!</div>
                    )}

                  {/* </Form.Item>
                </Form.Item> */}
              </div>
            )
          }
          <div className="chose">
            <Radio.Group onChange={onChangeRadio} defaultValue="sequence" buttonStyle="solid">
              <Radio.Button value="sequence">Sequence</Radio.Button>
              <Radio.Button value="fasta">Fasta</Radio.Button>
            </Radio.Group>
          </div>
          <div style={{display: 'flex', justifyContent:'center', alignItems: 'center'}}>
            {
              // result && {
                result?.map(item => (
                  <div
                    key={item.ID}
                    style={{color: '#fff', fontSize: '20px', width: '200px', margin: 'auto', border: 'solid 1px #ddd', borderRadius: '10px', marginBottom: '10px'}}
                  >
                    <div><span style={{fontWeight: 'bold'}}>Class:</span> {item.Class}</div>
                    <div><span style={{fontWeight: 'bold'}}>ID:</span> {item.ID}</div>
                    <div><span style={{fontWeight: 'bold'}}>Seq:</span> {item.Seq}</div>
                  </div>))
              // }
            }
          </div>
          <Button
            type="primary"
            icon={!loading && <CheckCircleOutlined />}
            htmlType="submit"
            style={{fontSize: '20px', height: 'auto', margin: '20px', width: '200px'}}
            disabled={loading}
          >
            {loading ? (<LoadingOutlined />) : ('Submit')}
          </Button>
        </Form>
      </div>
    </div>
  );
}

export default App;

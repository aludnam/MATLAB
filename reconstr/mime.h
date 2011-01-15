#ifndef MIME_H
#define MIME_H
#include <iostream>
#include <string>  // STL string class

using namespace std;


string extractTag(string & line, char * tag) // extracts a tag in format tag="value"
{
  string tags(tag);
  int offset=line.find(tag);
  int offset2=line.find("\"",offset+tags.length());
  string subs(line,offset+tags.length(),offset2-offset-tags.length());
  return subs;
}

void logline(const string & line)
{
  static ofstream *of=0;
  if (!of) of=new ofstream("upload/logfile.raw");
  if (of) 
    (*of) << line << flush;
  if (line.compare("closefile") == 0)
      of->close();
}
 
string ReadMime(istream& input, const string mimetag, string & content, string & filename)  // returns tagname value in name="value"
{
  string line2,empty,value;
  content=value;  // clear content
  filename=value;  // clear filename
  getline(input,line2);
  // logline(line2);
  // cout << mimetag << "\n"; 
  // cout << line2 << "\n"; 
  getline(input,empty);  // read empty line
  // logline(empty);
  getline(input,value);
  // logline(value);
  value += "\n"; // put the newline in again !
  // fill the content from the input until next mimetag
  while (value.compare(0,mimetag.length()-1,mimetag,0,mimetag.length()-1) != 0) {
    // logline("VALUE  :"+value+"\n");
    // logline("MIMETAG:"+mimetag+"\n");
    content+= value; // add string to content
    getline(input,value);
    value += "\n"; // put the newline in again !
    if (input.eof()) return "eof";
  }
  filename=extractTag(line2,"filename=\"");
  // logline("Filename:"+filename+"\n");

  return extractTag(line2,"name=\"");
}

#endif

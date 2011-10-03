/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef __TAB_FILE_WRITER_H__
#define __TAB_FILE_WRITER_H__

#include "MSReadRecord.h"
#include "TextFileWriter.h"

class TabFileWriter : public TextFileWriter {
 public:
  TabFileWriter(const std::string& filename);
  virtual ~TabFileWriter();
  void WriteRecord(const MSReadRecord& msread);
};

#endif /* __TAB_FILE_WRITER_H__ */

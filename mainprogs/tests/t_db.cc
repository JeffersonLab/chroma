// $Id: t_db.cc,v 1.1 2008-07-24 02:51:13 edwards Exp $
/*! \file
 *  \brief Test the Wilson mesons() routine
 */

#include "chroma.h"
#include "ConfDataStoreDB.h"

namespace Chroma
{
  using namespace FFDB;

  // Simple concrete key class
  class TestDBKey : public DBKey
  {
  public:
    // Part of Serializable
    const unsigned short serialID (void) const {return 789;}
    void objectBuffer (char* &buffer, unsigned int& length) {
      buffer = (char*)&data_; length = 2;
    }

    const unsigned int sizeOfObject (void) const throw (SerializeException) {
      sizeof(data_);
    }

    void writeObject (void) throw (SerializeException) {
      data_[0] = 1;
      data_[1] = 17;
    }

    void readObject (void) throw (SerializeException) {
    }

    // Part of DBKey
    int hasHashFunc (void) const {return 0;}
    int hasCompareFunc (void) const {return 0;}

    /**
     * Static Hash Function Implementation
     */
    static unsigned int hash (Db *db, const void* bytes, unsigned int len) {return 0;}

    /**
     * Static empty compare function 
     */
    static int compare (Db *db, const Dbt* k1, const Dbt* k2) {return 0;}

  private:
    int a;
    int b;
    int data_[2];
  };


  // Simple concrete data class
  class TestDBData : public DBData
  {
  public:
    const unsigned short serialID (void) const {return 789;}
    void objectBuffer (char* &buffer, unsigned int& length) {
      buffer = (char*)&data_; length = sizeOfObject();
    }

    const unsigned int sizeOfObject (void) const throw (SerializeException) {
      sizeof(data_);
    }

    void writeObject (void) throw (SerializeException) {
      for(int i=0; i < 10; ++i)
	data_[i] =  i;
    }

    void readObject (void) throw (SerializeException) {
    }


  private:
    float data_[10];
  };

}


using namespace Chroma;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_db.xml");
  push(xml,"t_db");

  push(xml,"lattice");
  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"Ns", Ns);
  write(xml,"nrow", nrow);
  pop(xml);

  // Try out some simple DB stuff
  {
    TestDBKey testDBKey;
    TestDBData testDBData;

    ConfDataStoreDB<TestDBKey, TestDBData> db("test.db");
  }

  // Done
  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}


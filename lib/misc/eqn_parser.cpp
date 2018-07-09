/**
 *  \file eqn_parser.cpp
 *
 *  \brief A string parser for equations
 *
 *  There is no precedence implied in here so parenthesis
 *  need to be used appropriatly, otherwise the evaluationis
 *  done left to right
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id:$
 */

/* System includes */
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctype.h>
#include<stdexcept>

#define DEF_MINVAL (0.0)
#define DEF_MAXVAL (1e35)

/* User includes */
extern "C" {
#include "ems.h"
#include "eqn_parser.h"
}

/**
 * Enum to denote various operations
 */
typedef enum {
  OP_NONE = -1,   /**< no operator*/
  OP_DONE,        /**< done, i.e. end of string */  
  OP_PLUS,        /**< operator+ */  
  OP_MINUS,       /**< operator- */  
  OP_MULT,        /**< operator* */  
  OP_DIV,         /**< operator/ */  
  OP_POW,         /**< operator^ */  
  OP_EXP,         /**< operator(exponential) */  
  OP_LOG,         /**< operator(natural logarithm) */  
  OP_NEG          /**< operator(neagte) */  
} operatorType;


/**
 * \brief An abstract base class for almost all classes in this parser
 *
 * This is the simplest block
 */
class eqnParserNode {
public:
  /**
   * Default destructor
   */
  virtual ~eqnParserNode(void) { };

  /**
   * The most important method. This needs to be virtual so the
   * most appropriate method is always called.
   * 
   * Performs the actual calcualation
   *
   * @return double value
   */
  virtual double getValue(void) = 0;

  /**
   * Gets the stored display string
   *
   * @param istr string to copy display string into (must be allocated)
   */
  virtual void getDisplayStr(char *istr) {
    strcpy(istr, dispStr);
  }

  /**
   * Sets (copies) the display string
   *
   * @param istr input display string, gets copied
   */
  void setDisplayStr(const char *istr) {
    strcpy(&dispStr[0], istr);
  }

private:
  /**
   * Holds the display string for this node
   */
  char dispStr[MAXSTRLEN];
};


/**
 * \brief Node to hold double values. i.e. constants in the equation
 */
class eqnParserNodeValue : public eqnParserNode {
public:
  /**
   * Default initialising constructor 
   *
   * @param val double value to store
   */
  eqnParserNodeValue(double val) {
    value = val;
  }

  /**
   * Get method for the value
   *
   * @return The double value
   */
  inline double getValue(void) {
    return(value);
  }

private:
  /**
   * Stored value
   */
  double value;
};


/**
 * \brief Node to hold pointers to value. i.e. variables
 */
class eqnParserNodeValuePtr : public eqnParserNode {
public:
  /**
   * Default constructor 
   *
   * @param valPtr double pointer to store
   */
  eqnParserNodeValuePtr(double *valPtr) {
    valuePtr = valPtr;
  }

  /**
   * Gets the value by dereferencing the pointer
   *
   * @return The double value
   */
  inline double getValue(void) {
    return(*valuePtr);
  }

private:
  /**
   * Stored pointer
   */
  const double *valuePtr;
};


/**
 * \brief Base Class to hold operations (unary)
 *
 * This class is what everything is cast'd to. Defines the operation
 * as well as the node (for the value) to operate on
 */
class eqnParserOp : public eqnParserNode {
public:
  /**
   *  Default constructor 
   */
  eqnParserOp(void) {
    setNode(NULL);
    minVal = DEF_MINVAL;
    maxVal = DEF_MAXVAL;
  }

  /**
   * Initialising constructor 
   *
   * @param inode node to assign to this operator instance
   */
  eqnParserOp(eqnParserNode *inode) {
    setNode(inode);
  }

  /* Destructor */
  ~eqnParserOp(void) {
    if (node != NULL)
      delete(node);
  }

  /* Set error string */
  void setErrStr(char *str) {
    strcpy(&errStr[0], str);
  }

  /* Get error string */
  char *getErrStr(void) {
    return(errStr);
  }

  /* Min/Max specification */
  void setMinMax(char **flds)
  {
    char *tag = flds[0];
    char *str = flds[1];
    double val;

    /* Set min value */
    if (strstr(tag, "min") && str2double(str, &val))
      minVal = val;
    
    /* Set max value */
    if (strstr(tag, "max") && str2double(str, &val))
      maxVal = val;
  }

  /** minVal getter
   *
   */
  inline double getMinVal(void) {
    return(minVal);
  }

  /** maxVal getter
   *
   */
  inline double getMaxVal(void) {
    return(maxVal);
  }

  /**
   * node getter
   *
   */
  inline eqnParserNode *getNode(void) {
    return(node);
  }

  /**
   * node setter
   */
  inline void setNode(eqnParserNode *inode) {
    node = inode;
  }
  
  /**
   * getter for the node value
   */
  inline double getNodeValue(void) {
    return(node->getValue());
  }

  /**
   * This class only has one node, so return the value of its node
   */
  virtual inline double getValue(void) {
    return(getNodeValue());
  }

  /**
   * Returns the operator string
   * NULL to indicate simple node
   */
  virtual const char *getOpStr(void) {
    return(NULL);
  }

  /**
   * Builds up the display string for this instance. Wraps parenthesis
   * for unary operators
   *
   * @param istr input [allocated] string
   */
  virtual void getDisplayStr(char *istr) {
    const char *opStr = getOpStr();
    /*
     * See there is a unary operator 
     */
    if (opStr != NULL) {
      sprintf(istr, "%s(", opStr);
      istr += strlen(opStr) + 1;
      getNode()->getDisplayStr(istr);
      istr += strlen(istr);
      sprintf(istr, ")");
    } else {
      /* Just copy from the node */
      getNode()->getDisplayStr(istr);
    }
  }

  /**
   * The display string to the left of this node. In this cass, there
   * are no neighbours, returnn itself
   */
  virtual void getLeftDisplayStr(char *istr) {
    getDisplayStr(istr);
  }

  /**
   * Actually print to stdout
   */
  char *display(void) {
    static char dispStr[MAXSTRLEN];
    getDisplayStr(&dispStr[0]);
    
    /* std::cout << dispStr << std::endl;*/
    return(dispStr);
  }

private:
  /**
   * Associated node object
   */
  eqnParserNode *node;

  /**
   * Error string
   */
  char errStr[MAXSTRLEN];

  /**
   * Min/max
   */
  double minVal;
  double maxVal;
};


/**
 * \brief Base class for all binary operations
 */
class eqnParserBiOp : public eqnParserOp {
public:
  /**
   * destructor 
   */
  ~eqnParserBiOp(void) {
    if (left_link != NULL)
      delete(left_link);
  }

  /**
   * Gets the value of the left link 
   */
  inline double getLeftLinkValue(void) {
    return(left_link->getValue());
  }

  /**
   * setter for the left link
   */
  inline void setLeftLink(eqnParserOp *ilink) {
    left_link = ilink;
  }

  /**
   * Builds up the display string from the left link leftwards
   */
  void getLeftDisplayStr(char *istr) {
    char buf[MAXSTRLEN];

    /* Get left value */
    left_link->getLeftDisplayStr(buf);
      
    /* Print the left value */
    sprintf(istr, "%s", buf);
    istr += strlen(buf);
      
    /* Add the operator */
    sprintf(istr, "%s", getOpStr());
    istr += strlen(getOpStr());

    /* Now get and print this value */
    getNode()->getDisplayStr(buf);
    sprintf(istr, "%s", buf);
    istr += strlen(buf);
  }

  /**
   * Builds up the display string for this node, wraps parenthesis
   */
  void getDisplayStr(char *istr) {
    char buf[MAXSTRLEN];

    sprintf(istr++, "(");

    /* Get left value */
    left_link->getLeftDisplayStr(buf);
      
    /* Print the left value */
    sprintf(istr, "%s", buf);
    istr += strlen(buf);
      
    /* Add the operator */
    sprintf(istr, "%s", getOpStr());
    istr += strlen(getOpStr());

    /* Now get and print this value */
    getNode()->getDisplayStr(buf);
    sprintf(istr, "%s", buf);
    istr += strlen(buf);

    /* Close out */
    sprintf(istr, ")");
  }

private:
  /**
   * Link to the left node
   */
  eqnParserOp *left_link;
};


/***************************/
/* Various implementations */
/***************************/
/**********/
/* BINARY */
/**********/

/**
 * \brief Addition
 */
class eqnParserPlus : public eqnParserBiOp {
public:
  /* Default constructor */
  eqnParserPlus(eqnParserNode *inode) {
    setNode(inode);
  }
  inline double getValue(void) {
    return(getLeftLinkValue() + getNodeValue());
  }
  inline const char * getOpStr(void) {
    return("+");
  }
};

/**
 * \brief Subtraction
 */
class eqnParserMinus : public eqnParserBiOp {
public:
  /* Default constructor */
  eqnParserMinus(eqnParserNode *inode) {
    setNode(inode);
  }
  inline double getValue(void) {
    return(getLeftLinkValue() - getNodeValue());
  }
  inline const char * getOpStr(void) {
    return("-");
  }
};

/**
 * \brief Multiplication
 */
class eqnParserMult : public eqnParserBiOp {
public:
  /* Default constructor */
  eqnParserMult(eqnParserNode *inode) {
    setNode(inode);
  }
  inline double getValue(void) {
    return(getLeftLinkValue() * getNodeValue());
  }
  inline const char * getOpStr(void) {
    return("*");
  }
};


/**
 * \brief Division
 */
class eqnParserDiv : public eqnParserBiOp {
public:
  /* Default constructor */
  eqnParserDiv(eqnParserNode *inode) {
    setNode(inode);
  }
  inline double getValue(void) {
    double rval = getNodeValue();
    if (rval != 0.0)
      return(getLeftLinkValue() / rval);
    else
      throw std::runtime_error("Divide by zero");
  }
  inline const char * getOpStr(void) {
    return("/");
  }
};


/**
 * \brief Power
 */
class eqnParserPow : public eqnParserBiOp {
public:
  /* Default constructor */
  eqnParserPow(eqnParserNode *inode) {
    setNode(inode);
  }
  inline double getValue(void) {
    return(pow(getLeftLinkValue(), getNodeValue()));
  }
  inline const char * getOpStr(void) {
    return("^");
  }
};

/*********/
/* UNARY */
/*********/
/**
 * \brief Exponential
 */
class eqnParserExp : public eqnParserOp {
public:
  /* Default constructor */
  eqnParserExp(eqnParserNode *inode) {
    setNode(inode);
  }
  inline double getValue(void) {
    return(exp(getNodeValue()));
  }
  inline const char * getOpStr(void) {
    return("exp");
  }
};


/**
 * \brief Logarithm
 */
class eqnParserLog : public eqnParserOp {
public:
  /* Default constructor */
  eqnParserLog(eqnParserNode *inode) {
    setNode(inode);
  }
  inline double getValue(void) {
    double ival = getNodeValue();
    if (ival > 0.0)
      /* Positive number, all good */
      return(log(ival));
    else {
      if (ival < 0)
	throw std::runtime_error("Log of negative number");
      else
	throw std::runtime_error("Log of zero");
    }
  }
  inline const char * getOpStr(void) {
    return("log");
  }
};


/**
 * \brief Negate
 */
class eqnParserNeg : public eqnParserOp {
public:
  /* Default constructor */
  eqnParserNeg(eqnParserNode *inode) {
    setNode(inode);
  }
  inline double getValue(void) {
    return(-(getNodeValue()));
  }
  inline const char * getOpStr(void) {
    return("-");
  }
};

/**
 * \brief Entry class used for the parsing
 */
class eqnParser {
public:
  /**
   * Default constructor 
   */
  eqnParser(void) {
    rem_str    = NULL;
    sub_str[0] = '\0';
    
    /**
     * Function pointers to handle custom variables
     */
    custom_data          = NULL;
    CustomGetValuePtrFcn = NULL;
  }
  
  /**
   * Initialising constructor 
   * @param istr String to parse
   */
  eqnParser(const char *istr) {
    rem_str    = NULL;
    sub_str[0] = '\0';

    /**
     * Function pointers to handle custom variables
     */
    custom_data          = NULL;
    CustomGetValuePtrFcn = NULL;

    setRemStr(istr);
  }
  
  
  /**
   * Destructor 
   */
  ~eqnParser(void) {
    delete[] orig_str;
  }

  /** 
   * getter for the original string that was passed in to parse
   */
  const char *getOrigStr(void) {
    return(orig_str);
  }

  /**
   * setter for the remaining string
   */
  inline void setRemStr(const char *istr) {
    rem_str = istr;
    /* Keep a copy for error reporting */
    orig_str = new char[strlen(rem_str)+1];
    strcpy(orig_str, istr);
  }

  /**
   * setter for the custom function
   */
  inline void setCustomFcn(eqn_custom_fcn fcn, void *data) {
    custom_data          = data;
    CustomGetValuePtrFcn = fcn;
  }
  
  /**
   * Marches string forward until the first non-blank character
   */
  void marchRemStr(int n) {
    int i=0;
    while(i++ < n)
      rem_str++;
  }
  
  /**
   * Create nodes based on whether they are a value, valuePtr or unary
   * operator 
   */
  eqnParserNode *createLeafNode(const char *tok) {
    double         val     = -1;
    double        *valPtr  = NULL;
    eqnParserNode *node    = NULL;
    
    /* See if its a double */
    if (str2double((char *)tok, &val)) {
      /* Yep ! */
      node = new eqnParserNodeValue(val);
    } else if ( (valPtr = CustomGetValuePtrFcn(tok, custom_data)) != NULL ) {
      node = new eqnParserNodeValuePtr(valPtr);
    } else {
      if (strlen(tok) > 1) {
	quit((char*)"(lib:misc:eqn_parser) Invalid variable name '%s' found\n", tok);
      } else {
	quit((char*)"(lib:misc:eqn_parser) Malformed equation found\n", 
	     (char *)tok);
      }
    }

    node->setDisplayStr(tok);

    return(node);
  }

  /** 
   * Parses the string until the next token and sets up the current
   * states
   */
  eqnParserOp *equate(void) {
    eqnParserOp *curr_lnk = NULL;
    eqnParserOp *prev_lnk = NULL;
	
    operatorType curr_op = OP_NONE;
    operatorType next_op = OP_NONE;
    operatorType un_op   = OP_NONE;

    int  curr_i = 0;
    char curr_tok[MAXSTRLEN];
    curr_tok[0] = '\0';

    /* Loop until the end of string */
    while (next_op != OP_DONE) {
      eqnParserNode *val_node = NULL;
      
      /* Skip over any spaces */
      while (isspace(*rem_str)) {
	rem_str++;
      }

      /*
       * Check if we've hit a unary operator
       */
      if ( (un_op = parseUnOpType(curr_tok)) > OP_NONE) {
	eqnParserNode *rec_node = NULL;
	eqnParserOp   *un_node  = NULL;

	/* Next character must be an open paren */
	if (*rem_str != '(') {
	  quit((char*)"(lib:misc:eqn_parser) unary operators must be immediatly followed by an open paren '%s'\n", curr_tok);
	}
	
	/* Recurse to create the underlying node */
	rec_node = doRecursion();
	
	/*
	 * Switch on the binary operator type create the OpNode
	 */
	switch(un_op) {
	case OP_EXP:
	  un_node = new eqnParserExp(rec_node);
	  break;
	case OP_LOG:
	  un_node = new eqnParserLog(rec_node);
	  break;
	case OP_NEG:
	  un_node = new eqnParserNeg(rec_node);
	  break;
	default:
	  quit((char*)"(lib:misc:eqn_parser) Unknown unary operator '%s' found\n", curr_tok);
	  break;
	}
	
	/* Reset curr_tok */
	curr_i = 0;
	curr_tok[0] = '\0';

	/* Make this the current value node */
	val_node = un_node;
      }
      
      /*
       * Recurse at parenthesis
       */
      if (*rem_str == '(') {
	val_node = doRecursion();
      }

      /* See if we've hit a binary operator or reached the end */
      if ( ((next_op = parseBiOpType()) > OP_NONE) ) {
	/*
	 * Create this leaf node, could be unary
	 */
	if (val_node == NULL)
	  val_node = createLeafNode(curr_tok);
	
	/* See if this is the very frist node */
	if (curr_op == OP_NONE) {
	  /* This is the very first node */
	  curr_lnk = new eqnParserOp(val_node);

	} else {
	  eqnParserBiOp *curr_bi_lnk = NULL;

	  /*
	   * Switch on the binary operator type create the OpNode
	   */
	  switch (curr_op) {
	  case OP_PLUS:
	    curr_bi_lnk = new eqnParserPlus(val_node);
	    break;
	  case OP_MINUS:
	    curr_bi_lnk = new eqnParserMinus(val_node);
	    break;
	  case OP_MULT:
	    curr_bi_lnk = new eqnParserMult(val_node);
	    break;
	  case OP_DIV:
	    curr_bi_lnk = new eqnParserDiv(val_node);
	    break;
	  case OP_POW:
	    curr_bi_lnk = new eqnParserPow(val_node);
	    break;
	  default:
	    quit((char*)"(lib:misc:eqn_parser) Unknown binary operator '%s' found\n", curr_tok);
	    break;
	  }

	  /* Chain them up */
	  if (prev_lnk != NULL) {
	    curr_bi_lnk->setLeftLink(prev_lnk);
	  }
	  curr_lnk = (eqnParserOp *)curr_bi_lnk;
	}

	/* Reset current token counter and op */
	curr_i      = 0;
	curr_tok[0] = '\0';
	curr_op     = next_op;

	prev_lnk = (eqnParserOp *)curr_lnk;
	
      } else {
	/* Keep caching the current token */
	curr_tok[curr_i++] = *rem_str;
	curr_tok[curr_i]   = '\0';
      }

      /*
       * Move forward one.  I'm not convinced that this is the best
       * place for this increment as it should reside within the type of
       * binary operator but I haven't been able to nail down the most
       * elegant solution yet.
       */
      marchRemStr(1);
    }

    /* Return the last node */
    return(curr_lnk);
  }

  
  /**
   * Creates the substring with parenthesis
   */
  eqnParserNode *doRecursion(void) {
    int str_i = 0;
    int stack = 1; /* We're already within a paren */
    bool done = false;
    eqnParserOp *node = NULL;

    /* Create object */
    eqnParser eqn;

    /* Go past the open paren */
    rem_str++;

    /* March along and recurse, if needed */
    while ( (!done && *rem_str != '\0') ) {

      /* Match up the parens */
      if (*rem_str == '(') {
	/* push */
	stack++;
      } else if (*rem_str == ')') {
	/* pop */
	stack--;
      }

      /* Check for empty stack */
      if (stack == 0)
	done = true;
      else
	sub_str[str_i++] = *rem_str;

      /* Move Along */
      rem_str++;

      /* Skip over any spaces */
      while (isspace(*rem_str)) {
	rem_str++;
      }
    }
  
    /* Terminate the string */
    sub_str[str_i] = '\0';
    
    /* Set the string and ask for the child node */
    eqn.setRemStr(&sub_str[0]);

    /* Set up the custom data */
    eqn.setCustomFcn(CustomGetValuePtrFcn, custom_data);

    /* Parses and returns the last node */
    node = eqn.equate();

    return(node);
  }

  /**
   * Returns operatorType based on the binary operator character
   */
  inline operatorType parseBiOpType(void) {
    operatorType op;

    if (*rem_str == '+') {
      op = OP_PLUS;
    } else if (*rem_str == '-') {
      op = OP_MINUS;
    } else if (*rem_str == '*') {
      op = OP_MULT;
    } else if (*rem_str == '/') {
      op = OP_DIV;
    } else if (*rem_str == '^') {
      op = OP_POW;
    } else if (*rem_str == '\0') {
      op = OP_DONE;
    } else
      op = OP_NONE;

    return(op);
  }

  /**
   * Returns operatorType based on the unary operator keyword
   */
  inline operatorType parseUnOpType(const char *itok) {
    operatorType op;

    if (strcmp(itok, "exp") == 0) {
      op = OP_EXP;
    } else if (strcmp(itok, "log") == 0) {
      op = OP_LOG;
    } else if (strcmp(itok, "neg") == 0) {
      op = OP_NEG;
    } else
      op = OP_NONE;
    
    return(op);
  }

private:
  /**
   * This is the remaining string 
   */
  const char *rem_str;

  /**
   * This is a buffer for temp substrings 
   */
  char sub_str[MAXSTRLEN];

  /**
   * Original string for error reporting 
   */
  char *orig_str;

  /**
   * Custom data pointer provided by the user to be used as a callback
   * for getting the value pointers
   *
   * Note: memory management of this pointer *must* be handled by the caller
   */
  void *custom_data;
  eqn_custom_fcn CustomGetValuePtrFcn;
};

/******************/
/* Public C api's */
/******************/ 

/**
 * Main parser constructor
 *
 * @param str the string to parse
 * @param fcn function for the parser to call for unknown tokens
 * @param data custom data for use with fcn
 * @param err Error string to use for division by zero errors
 */
extern "C"
void *EqnCreateParser(const char *str, eqn_custom_fcn fcn, void *data, char *err)
{
  int nf;
  char *fields[MAXSTRLEN * 3];
  char *lstr = strdup(str);

  /* Create main object */
  eqnParser node;
  
  /* This is the last node */
  eqnParserOp *p = NULL;

  /* See how many semi-colon tokens we have */
  nf = split(lstr, fields, (char *)";");
  
  /* Set the last field as remaining string */
  node.setRemStr(fields[nf-1]);

  /* Set the custom callback */
  node.setCustomFcn(fcn, data);

  /* Do the parsing */
  p = node.equate();

  /* Check for min/max specification */
  while (--nf) {
    char *f2[MAXSTRLEN * 2];
    int nf2 = split(fields[nf-1], f2, (char *)"=");
    if (nf2 != 2)
      quit((char *)"(lib:misc:eqn_parser) Malformed min/max specification '%s'\n", str);
    p->setMinMax(f2);
  }

  /* Set the error string */
  p->setErrStr(err);

  /* Return as a void pointer */
  return((void *)p);
}

/**
 * Returns the display string
 *
 * @param e pointer to the parser object created by EqnCreateParser
 * @return the string
 */
extern "C"
char *EqnDisplayStr(void *e)
{
  double min_val = ((eqnParserOp *)e)->getMinVal();
  double max_val = ((eqnParserOp *)e)->getMaxVal();
  
  char *str = ((eqnParserOp *)e)->display();

  /* Append any bounds */
  if (min_val != DEF_MINVAL)
    sprintf(&str[strlen(str)], "; min = %f", min_val);

  if (max_val != DEF_MAXVAL)
    sprintf(&str[strlen(str)], "; max = %f", max_val);

  return(str);
}

/**
 * Evaluates the value of this equation
 * Note: Clamps negative values to zero
 *
 * @param e pointer to the parser object created by EqnCreateParser
 * @return the double value
 */
extern "C"
double EqnGetValue(void *e)
{
  double min_val = ((eqnParserOp *)e)->getMinVal();
  double max_val = ((eqnParserOp *)e)->getMaxVal();

  try {
    /* Get the value, may throw */
    double val = ((eqnParserOp *)e)->getValue();

    /* Return within bounds */
    return(val < min_val ? min_val : (val > max_val ? max_val : val));

  } catch (std::runtime_error &re) {
    /* Catch run time errors */
    quit((char *)"(lib:misc:eqn_parser) %s when evaluating %s!\n",
	 re.what(), ((eqnParserOp *)e)->getErrStr());
  }
  
  return(0);
}

/**
 * De-allocate memory associated with the parser object
 *
 * @param e pointer to the parser object created by EqnCreateParser
 */
extern "C"
void EqnFree(void *e)
{
  if (e != NULL)
    delete((eqnParserOp *)e);
}

/* EOF */

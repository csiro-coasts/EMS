
/* 
------------------------------------------------------------
 * DEL-TRI-TYPE-2: data-types for 2-dim. delaunay tria.
------------------------------------------------------------
 *
 * This program may be freely redistributed under the 
 * condition that the copyright notices (including this 
 * entire header) are not removed, and no compensation 
 * is received through use of the software.  Private, 
 * research, and institutional use is free.  You may 
 * distribute modified versions of this code UNDER THE 
 * CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
 * TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
 * ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
 * MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
 * NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
 * of this code as part of a commercial system is 
 * permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
 * AUTHOR.  (If you are not directly supplying this 
 * code to a customer, and you are instead telling them 
 * how they can obtain it for free, then you are not 
 * required to make any arrangement with me.) 
 *
 * Disclaimer:  Neither I nor: Columbia University, The
 * Massachusetts Institute of Technology, The 
 * University of Sydney, nor The National Aeronautics
 * and Space Administration warrant this code in any 
 * way whatsoever.  This code is provided "as-is" to be 
 * used at your own risk.
 *
------------------------------------------------------------
 *
 * Last updated: 12 May, 2017
 *
 * Copyright 2013-2017
 * Darren Engwirda
 * de2363@columbia.edu
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __DELAUNAY_TRI_TYPE_2__
#   define __DELAUNAY_TRI_TYPE_2__

    namespace mesh {

/*
------------------------------------------------------------
 * DELAUNAY-TRI-NODE-2: node for delaunay tria. in R^2.
------------------------------------------------------------
 * IPTR-TYPE - signed-integer typedef.
 * REAL-TYPE - floating-point typedef.
------------------------------------------------------------
 */
     
    template <
    typename I,
    typename R
             >
    class delaunay_tri_node_2
    {
/*----------- node-type for delaunay triangulation in R^2 */
    public  :
    typedef R               real_type ;
    typedef I               iptr_type ;

    iptr_type static constexpr  _dims = +2 ;

/*------------------------------------ local data members */
    containers::
    fixed_array<real_type, +2>  _pval ;  // node coord
    
    containers::
    fixed_array<char_type, +4>  _flag ;  // kind, mark
    
    iptr_type                   _next ;

/*---------------------------------------- "write access" */
    __inline_call real_type &pval (
        iptr_type  _ipos
        )
    {   return this->_pval [_ipos];
    }
    __inline_call iptr_type &next (
        )
    {   return this->_next;
    }
    __inline_call char_type &mark (
        )
    {   return this->_flag [   +0];
    }

/*---------------------------------------- "const" access */
    __inline_call real_type const& pval (
        iptr_type _ipos
        ) const
    {   return this->_pval [_ipos];
    }
    __inline_call iptr_type const& next (
        ) const
    {   return this->_next;
    }
    __inline_call char_type const& mark (
        ) const
    {   return this->_flag [   +0];
    }

    } ;

/*
------------------------------------------------------------
 * DELAUNAY-TRI-TRIA-2: tria for delaunay tria. in R^2.
------------------------------------------------------------
 * IPTR-TYPE - signed-integer typedef.
 * REAL-TYPE - floating-point typedef.
------------------------------------------------------------
 */
 
    template <
    typename I,
    typename R
             >
    class delaunay_tri_tria_2
    {
/*----------- tria-type for delaunay triangulation in R^2 */
    public  :
    typedef R               real_type ;
    typedef I               iptr_type ;

    iptr_type static constexpr  _dims = +2 ;

/*------------------------------------ local data members */
    containers::
    fixed_array<iptr_type, +3>  _ndat ;  // node indexing
    
    containers::
    fixed_array<iptr_type, +3>  _edat ;  // edge neighbour
    
    containers::
    fixed_array<char_type, +3>  _epos ;  // edge neighbour
    
    containers::
    fixed_array<char_type, +4>  _flag ;  // kind, mark

/*---------------------------------------- "write" access */
    __inline_call iptr_type &node (
        iptr_type  _ipos
        )
    {   return this->_ndat [_ipos];
    }
    __inline_call iptr_type &next (
        iptr_type  _ipos
        )
    {   return this->_edat [_ipos];
    }
    __inline_call char_type &fpos (
        iptr_type  _ipos
        )
    {   return this->_epos [_ipos];
    }
    __inline_call char_type &mark (
        )
    {   return this->_flag [   +0];
    }

/*---------------------------------------- "const" access */
    __inline_call iptr_type const& node (
        iptr_type  _ipos
        ) const
    {   return this->_ndat[_ipos];
    }
    __inline_call iptr_type const& next (
        iptr_type  _ipos
        ) const
    {   return this->_edat[_ipos];
    }
    __inline_call char_type const& fpos (
        iptr_type  _ipos
        ) const
    {   return this->_epos[_ipos];
    }
    __inline_call char_type const& mark (
        ) const
    {   return this->_flag[   +0];
    }

/*---------------------------------------- local indexing */
    __static_call
    __inline_call void_type faceind22 (
        iptr_type *_fnod,
        iptr_type//_fpos
        )
    {
        _fnod[0] = 0 ;
        _fnod[1] = 1 ;
        _fnod[2] = 2 ;
    }
    __static_call
    __inline_call void_type faceind11 (
        iptr_type *_fnod,
        iptr_type//_fpos
        )
    {
        _fnod[0] = 0 ;
        _fnod[1] = 1 ;
    }
    __static_call
    __inline_call void_type faceind21 (
        iptr_type *_fnod,
        iptr_type  _fpos
        )
    {
        switch (_fpos)
            {
        case 0 :
            {
        _fnod[0] = 0 ;
        _fnod[1] = 1 ; 
        _fnod[2] = 2 ; break ;
            }
        case 1 :
            {
        _fnod[0] = 1 ;
        _fnod[1] = 2 ; 
        _fnod[2] = 0 ; break ;
            }
        case 2 :
            {
        _fnod[0] = 2 ;
        _fnod[1] = 0 ; 
        _fnod[2] = 1 ; break ;
            }
            }
    }
    __static_call
    __inline_call void_type face_node  (
        iptr_type *_fnod,
        iptr_type  _fpos,
        iptr_type  _from,
        iptr_type  _into
        )
    {
    /* index FROM 2-dim faces */
        if (_from == +2)
            {
    /* index INTO 2-dim faces */
        if (_into == +2)
            {
        faceind22(_fnod, _fpos);
            }
    /* index INTO 1-dim faces */
        else
        if (_into == +1)
            {
        faceind21(_fnod, _fpos);
            }
            }
    /* index FROM 1-dim faces */
        else
        if (_from == +1)
            {
    /* index INTO 1-dim faces */
        if (_into == +1)
            {
        faceind11(_fnod, _fpos);
            }
            }
    }

    } ;


    }
    
#   endif//__DELAUNAY_TRI_TYPE_2__  
    
    
    

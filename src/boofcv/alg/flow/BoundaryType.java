package boofcv.alg.flow;

// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.
//
// Java port modifications (C) 2014 by Peter Abeles <peter.abeles@gmail.com>

/**
 * @author Peter Abeles
 */
public enum BoundaryType
{
   DIRICHLET,
   REFLECTING,
   PERIODIC;

   public static final BoundaryType DEFAULT = REFLECTING;
}

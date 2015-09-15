/* 
 * Lattice.java
 *
 * Defines a new "Lattice" type, which is a directed acyclic graph that
 * compactly represents a very large space of speech recognition hypotheses
 *
 * Note that the Lattice type is immutable: after the fields are initialized
 * in the constructor, they cannot be modified.
 *
 * Students may only use functionality provided in the packages
 *     java.lang
 *     java.util 
 *     java.io
 *     
 * as well as the class java.math.BigInteger
 * 
 * Use of any additional Java Class Library components is not permitted 
 * 
 * Nick Nestor    Nathan Hansen
 * W01147467      W01005100
 *
 */
 
 import java.lang.*;
 import java.util.*;
 import java.io.*;
 import java.math.BigInteger;
 
 
public class Lattice {
    private String utteranceID;       // A unique ID for the sentence
    private int startIdx, endIdx;     // Indices of the special start and end tokens
    private int numNodes, numEdges;   // The number of nodes and edges, respectively
    private Edge[][] adjMatrix;       // Adjacency matrix representing the lattice
                                      //   Two dimensional array of Edge objects
                                      //   adjMatrix[i][j] == null means no edge (i,j)
    private double[] nodeTimes;       // Stores the timestamp for each node



//-----------------------------------------LATTICE------------------------------------------//
// complete


    // Lattice
    // Preconditions:
    //     - latticeFilename contains the path of a valid lattice file
    // Post-conditions
    //     - Field id is set to the lattice's ID
    //     - Field startIdx contains the node number for the start node
    //     - Field endIdx contains the node number for the end node
    //     - Field numNodes contains the number of nodes in the lattice
    //     - Field numEdges contains the number of edges in the lattice
    //     - Field adjMatrix encodes the edges in the lattice:
    //        If an edge exists from node i to node j, adjMatrix[i][j] contains
    //        the address of an Edge object, which itself contains
    //           1) The edge's label (word)
    //           2) The edge's acoustic model score (amScore)
    //           3) The edge's language model score (lmScore)
    //        If no edge exists from node i to node j, adjMatrix[i][j] == null
    //
    //         IS THIS SUPPOSED TO BE HERE?
    //     - Field nodeTimes contains the timestamps of the nodes in the lattice
    //
    // Notes:
    //     - If you encounter a FileNotFoundException, print to standard error
    //         "Error: Unable to open file " + latticeFilename
    //       and exit with status (return code) 1
    //     - If you encounter a NoSuchElementException, print to standard error
    //         "Error: Not able to parse file " + latticeFilename
    //       and exit with status (return code) 2
    public Lattice(String latticeFilename) {
      Scanner scan = null;
      
      try {
         scan = new Scanner(new File(latticeFilename));
         }
      catch(FileNotFoundException ex) {
         System.out.println("Error: Unable to open file " + latticeFilename);
         System.exit(1);
         }
      
      // set the lattice ID to the second segment of the first line
      String id_set[] = ((String)scan.nextLine()).split("\\s");
      utteranceID = id_set[1];
      
      // set the starting index to the designated node in the lattice - line 2
      String sidx_set[] = ((String)scan.nextLine()).split("\\s");
      startIdx = Integer.parseInt(sidx_set[1]);
      
      // set the ending index to the designated node in the lattice - line 3
      String eidx_set[] = ((String)scan.nextLine()).split("\\s");
      endIdx = Integer.parseInt(eidx_set[1]);
      
      // set the number of nodes to the designated value in the lattice - line 4
      String nodeNum_set[] = ((String)scan.nextLine()).split("\\s");
      numNodes = Integer.parseInt(nodeNum_set[1]);
      
      // set the number of edges to the designated value in the lattice - line 5
      String edgeNum_set[] = ((String)scan.nextLine()).split("\\s");
      numEdges = Integer.parseInt(edgeNum_set[1]);
      
      // the array and matrix to hold all the node and edge information
      nodeTimes = new double[numNodes];
      adjMatrix = new Edge[numNodes][numNodes];
      
      while (scan.hasNext()) {
         // if next line describes a node
         if (scan.hasNext("node")) {
            String nodeHolder = (String)scan.next();
            int node = scan.nextInt();
            double time = scan.nextDouble();
            
            nodeTimes[node] = time;
            }
         // if next line describes an edge
         else if (scan.hasNext("edge")) {
            String edgeHolder = (String)scan.next();
            int i = scan.nextInt();
            int j = scan.nextInt();
            String label = (String)scan.next();
            int amScore = scan.nextInt();
            int lmScore = scan.nextInt();
            
            Edge edge = new Edge(label, amScore, lmScore);
            adjMatrix[i][j] = edge;
            }
         // the file has a formatting error. fix it.
         else {
            System.out.println("Error: file formatting issue, please follow the layout: ");
            System.out.print("'node nodeNumber nodeTimeStamp' or ");
            System.out.println("'edge sendingNode receivingNode edgeLabel amScore lmScore'");
            }
         }
      }   
    
    
//--------------------------------------GETUTTERANCEID---------------------------------------//    
// complete


    // getUtteranceID
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the utterance ID
    public String getUtteranceID() {
        return utteranceID;
      }
    
    
    
//--------------------------------------GETNUMNODES------------------------------------------//    
// complete


    // getNumNodes
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the number of nodes in the lattice
    public int getNumNodes() {
        return numNodes;
      }
    
    
    
//-------------------------------------GETNUMEDGES------------------------------------------//
// complete
    

    // getNumEdges
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the number of edges in the lattice
    public int getNumEdges() {
        return numEdges;
      }



//--------------------------------------TOSTRING------------------------------------------//
// complete


    // toString
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Constructs and returns a string that is identical to the contents
    //      of the lattice file used in the constructor
    // Notes:
    //    - Do not store the input string verbatim: reconstruct it on they fly
    //      from the class's fields
    //    - toString simply returns a string, it should not print anything itself
    // Hints:
    //    - You can use the format method on a PrintStream to print a floating
    //      point value to two decimal places
    public String toString() {
      StringBuilder returnString = new StringBuilder("id  " + utteranceID + "\n" +
         "start  " + startIdx + "\n" + "end " + endIdx + "\n" + "numNodes " + numNodes + "\n" +
         "numEdges " + numEdges + "\n");
      
      for (int node = 0; node < numNodes; node++) {
         returnString.append("node " + node + " " + nodeTimes[node] + "\n");
         }
         
      for (int sidx = 0; sidx < numNodes; sidx++) {
         for (int eidx = 1; eidx < numNodes; eidx++) {
            if (adjMatrix[sidx][eidx] != null) {
               returnString.append("edge " + sidx + " " + eidx + " "
               + (adjMatrix[sidx][eidx]).getLabel() + " " + (adjMatrix[sidx][eidx]).getAmScore()
               + " " + (adjMatrix[sidx][eidx]).getLmScore() + "\n");
               }
            }
         }
         
      return returnString.toString();
      }
    
    
    
//-----------------------------------------DECODE-------------------------------------------//    
// complete!


    // decode
    // Pre-conditions:
    //    - lmScale specifies how much lmScore should be weighted
    //        the overall weight for an edge is amScore + lmScale * lmScore
    // Post-conditions:
    //    - A new Hypothesis object is returned that contains the shortest path
    //      (aka most probable path) from the startIdx to the endIdx
    // Hints:
    //    - You can create a new empty Hypothesis object and then
    //      repeatedly call Hypothesis's addWord method to add the words and 
    //      weights, but this needs to be done in order (first to last word)
    //      Backtracking will give you words in reverse order.
    //    - java.lang.Double.POSITIVE_INFINITY represents positive infinity
    public Hypothesis decode(double lmScale) {
      Hypothesis shortestPath = new Hypothesis();
      double[] cost = new double[numNodes];
      int[] parent = new int[numNodes];
      int[] sortedNodes = topologicalSort();
      Stack<Integer> ordered = new Stack<Integer>(); 
      int node = endIdx;
      
      for (int j = 1; j < numNodes; j++) {
         cost[j] = java.lang.Double.POSITIVE_INFINITY;
         }
      
      cost[startIdx] = 0.0;
      
      for (int i = 0; i < numNodes; i++) {
         for (int j = 0; j < numNodes; j++) {
            if (adjMatrix[i][j] != null) {
               double weight = adjMatrix[i][j].getAmScore() + (lmScale * adjMatrix[i][j].getLmScore());
               
               // System.out.println("edge from " + i + " to " + j);
               
               if ((weight + cost[i]) < cost[j]) {
                  // System.out.println("edge from " + i + " to " + j);
                  // System.out.println("weight of " + (weight + cost[i]) + " is less than " + cost[j]);
                  cost[j] = (weight + cost[i]);
                  parent[j] = i;
                  }
               }
            }
         }
      
      while (node != startIdx) {
         ordered.push(node);
         node = parent[node];
      }
      
      ordered.push(node);
      
      while (!(ordered.isEmpty())) {
         node = ordered.pop();
         if (!(ordered.isEmpty())) {
            shortestPath.addWord(adjMatrix[node][ordered.peek()].getLabel(), cost[ordered.peek()]);
            }
         }
      
      return shortestPath;
      }
    
    
    
//----------------------------------TOPOLOGICALSORT------------------------------------------//
// complete

    
    // topologicalSort
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - A new int[] is returned with a topological sort of the nodes
    //      For example, the 0'th element of the returned array has no 
    //      incoming edges.  More generally, the node in the i'th element 
    //      has no incoming edges from nodes in the i+1'th or later elements
   public int[] topologicalSort() {
      int[] result = new int[numNodes];
      int resultIdx = 0;
      int[] inDegrees = new int[numNodes];
      int inDegreeSum = 0;
      HashSet<Integer> zeroInDegree = new HashSet<Integer>();
      HashSet<Integer> completed = new HashSet<Integer>();
      zeroInDegree.add(startIdx);
      
      // Initialize the inDegree array. Index corresponds to node number, and
      // value corresponds to in-degree.
      for (int i = 0; i < numNodes; i++) {
         for (int j = 1; j < numNodes; j++) {
            // System.out.println(adjMatrix[i][j]);
            if (adjMatrix[i][j] != null) {
               inDegrees[j] = inDegrees[j] + 1;
               // System.out.println("in-degree of " + j + " is " + inDegrees[j]);
               }
            }
         }
         
      // While there are still nodes with an in-degree of zero, find the lowest
      // node, add it to the result array, decrement the inDegree array to remove
      // the lowest node's edges, and add the recieving node if it now has 0 in-
      // degree.
      while (!zeroInDegree.isEmpty()) {
         int lowestNode = endIdx;
         // System.out.println("Setting lowestNode to endIdx: " + endIdx);
         Iterator<Integer> it = zeroInDegree.iterator();
         
         while (it.hasNext()) {
            int nextNode = it.next();
            // System.out.println("nextNode in zeroInDegree: " + nextNode);
            if (nextNode < lowestNode) {
               // System.out.println("nextNode [" + nextNode + "] is less than lowestNode ["+ lowestNode + "] ");
               lowestNode = nextNode;
               }
            }
         
         zeroInDegree.remove(lowestNode);
         // System.out.println("Removing node in zeroInDegree: " + lowestNode);
         // System.out.println("Setting node to result[" + resultIdx + "] ");
         result[resultIdx] = lowestNode;
         completed.add(lowestNode);
         resultIdx++;
         
         for (int j = 1; j < numNodes; j++) {
            if (adjMatrix[lowestNode][j] != null) {
               inDegrees[j] = inDegrees[j] - 1;
               // System.out.println("Decrementing in-degree of node " + j);
               }
            if (inDegrees[j] == 0 && !(zeroInDegree.contains(j)) && !(completed.contains(j))) {
               // System.out.println("Adding node to zeroInDegree: " + j);
               zeroInDegree.add(j);
               }
            }
         }
      
      // loop checking!
      for (int j = 1; j < numNodes; j++) {
         inDegreeSum += inDegrees[j];
         }
      
      if (inDegreeSum > 0) {
         System.out.println("Error: cycle found. ");
         }
      
      return result;
      }
         
      
      
      
      
    
    
    
//-----------------------------------COUNTALLPATHS------------------------------------------//
// INCORRECT
    

    // countAllPaths
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the total number of distinct paths from startIdx to endIdx
    // Hints:
    //    - The straightforward recursive traversal is prohibitively slow
    //    - This can be solved efficiently using something similar to the 
    //        shortest path algorithm used in decode
    //        Instead of min'ing scores over the incoming edges, you'll want to 
    //        do some other operation...
    public java.math.BigInteger countAllPaths() {
//       BigInteger numPaths = new BigInteger("0");
//       BigInteger incrementor = new BigInteger("1");
//       ArrayList<LinkedList<Integer>> adjSet = new ArrayList<LinkedList<Integer>>();
//       ArrayList<LinkedList<Integer>> visited = new ArrayList<LinkedList<Integer>>();
//       Stack<Integer> toVisit = new Stack<Integer>();     
//       
//       // create the adjSet that holds a node as index, and a linked list as the nodes it connects to
//       for (int i = 0; i < numNodes; i++) {
//          LinkedList<Integer> adjacentNodes = new LinkedList<Integer>();
//          
//          for (int j = 1; j < numNodes; j++) {
//             if (adjMatrix[i][j] != null) {
//                adjacentNodes.add(j);
//                }
//             }
//          
//          adjSet.add(adjacentNodes);
//          }
//       
//       toVisit.push(startIdx);
//       
//       // while there are still paths we should find
//       while (!toVisit.isEmpty()) {
//          int i = toVisit.peek();
//          int finished = 1;
//          int pushed = 0;
//          
//          // check if all of the nodes in the adjacency set have already been checked
//          for (int j : adjSet.get(i)) {
//             System.out.println("adjSet of [" + i + "] is " + adjSet.get(i));
//             System.out.println("Checking if node [" + j + "] is in the adjSet. ");
//             if (!visited.get(i).contains(j)) {
//                System.out.println("Setting finished to 0. \n");
//                finished = 0;
//                }
//             }
//          
//          // if there aren't any nodes left to visit, pop the node
//          if (finished == 1) {
//             toVisit.pop();
//             }
//          else {
//             // extract a single element from the adjacency set and visit it
//             for (int j : adjSet.get(i)) {
//                if (pushed != 1 && !(adjSet.get(i)).equals(visited.get(i))) {
//                   if (!visited.contains(j)) {
//                      toVisit.push(j);
//                      visited.get(i).add(j);
//                      }
//                   }
//                }
//             if (toVisit.peek() == endIdx) {
//                numPaths.add(incrementor);
//                visited.clear();
//                }
//             }
//          }
         
      return backwardsCount(endIdx);
      }
   
   private java.math.BigInteger backwardsCount (int j) {
      BigInteger paths = new BigInteger("0");
      BigInteger incrementor = new BigInteger("1");
      
      // for the sending nodes to this node, add all the possible paths
      for (int i = 0; i < numNodes; i++) {
         if (adjMatrix[i][j] != null) {
            paths = paths.add(backwardsCount(i));
            }
         }
      
      // no previous nodes, return 1 path
      if (j == startIdx) {
         return incrementor;
         }
      
      return paths;
      }
    
    
    
//--------------------------------GETLATTICEDENSITY------------------------------------------//
// complete
    

    // getLatticeDensity
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the lattice density, which is defined to be:
    //      (# of non -silence- words) / (# seconds from start to end index)
	 //      Note that multiwords (e.g. to_the) count as a single non-silence word
   public double getLatticeDensity() {
      double density = 0.0;
        
      for (int i = 0; i < numNodes; i++) {
         for (int j = 1; j < numNodes; j++) {
            if (!(adjMatrix[i][j] == null) && !(adjMatrix[i][j].getLabel().equals("-silence-"))) {
               density++;
               }
            }
         }
         
      return (density / (nodeTimes[endIdx] - nodeTimes[startIdx]));   
      }



//------------------------------------WRITEASDOT------------------------------------------//
// complete


    // writeAsDot - write lattice in dot format
    // Pre-conditions:
    //    - dotFilename is the name of the intended output file
    // Post-conditions:
    //    - The lattice is written in the specified dot format to dotFilename
    // Notes:
    //    - See the assignment description for the exact formatting to use
    //    - For context on the dot format, see    
    //        - http://en.wikipedia.org/wiki/DOT_%28graph_description_language%29
    //        - http://www.graphviz.org/pdf/dotguide.pdf
   public void writeAsDot(String dotFilename) {
      String arrow = " -> ";
      
      try {
			File dotFile = new File(dotFilename);
 
			if (!dotFile.exists())
				dotFile.createNewFile();
 
			FileWriter fileWriter = new FileWriter(dotFile.getAbsoluteFile());
			BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
         
         bufferedWriter.write("digraph g { \n\trankdir=\"LR\" \n");
         
         for (int i = 0; i < numNodes; i++) {
            for (int j = 1; j < numNodes; j++) {
               if (!(adjMatrix[i][j] == null)) {
                  bufferedWriter.write("\t" + i + arrow + j + " ");
                  bufferedWriter.write("[label = \"" + adjMatrix[i][j].getLabel() + "\"] \n");
                  }
               }
            }
            
         bufferedWriter.write("}");
            
         bufferedWriter.close();   
         } catch (IOException ex) {
			ex.printStackTrace();
		   }         
      }



//--------------------------------------SAVEASFILE------------------------------------------//
// complete


    // saveAsFile - write in the simplified lattice format (same as input format)
    // Pre-conditions:
    //    - latticeOutputFilename is the name of the intended output file
    // Post-conditions:
    //    - The lattice's toString() representation is written to the output file
    // Note:
    //    - This output file should be identical to the original input .lattice file
   public void saveAsFile(String latticeOutputFilename) {
      try {
			File outputFile = new File(latticeOutputFilename);
 
			if (!outputFile.exists())
				outputFile.createNewFile();
 
			FileWriter fileWriter = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
         
         String fileContents = toString();
         
         bufferedWriter.write(fileContents);
         bufferedWriter.close(); 
		   } catch (IOException ex) {
			ex.printStackTrace();
		   }
      }



//----------------------------------UNIQUEWORDSATTIME------------------------------------------//
// complete


    // uniqueWordsAtTime - find all words at a certain point in time
    // Pre-conditions:
    //    - time is the time you want to query
    // Post-conditions:
    //    - A HashSet is returned containing all unique words that overlap 
    //      with the specified time
    //     (If the time is not within the time range of the lattice, the Hashset should be empty)
   public java.util.HashSet<String> uniqueWordsAtTime(double time) { 
      HashSet<String> uniqueWords = new HashSet<String>();
        
      for (int i = 0; i < numNodes; i++) {
         for (int j = 1; j < numNodes; j++) {
            if (adjMatrix[i][j] != null && (nodeTimes[i] < time && nodeTimes[j] > time)) {
               uniqueWords.add(adjMatrix[i][j].getLabel());
               }
            }
         }
         
      return uniqueWords;
      }



//------------------------------------PRINTSORTEDHITS------------------------------------------//
// complete


    // printSortedHits - print in sorted order all times where a given token appears
    // Pre-conditions:
    //    - word is the word (or multiword) that you want to find in the lattice
    // Post-conditions:
    //    - The midpoint (halfway between start and end time) for each instance of word
    //      in the lattice is printed to two decimal places in sorted (ascending) order
    //      All times should be printed on the same line, separated by a single space character
    //      (If no instances appear, nothing is printed)
    // Note:
    //    - java.util.Arrays.sort can be used to sort
    //    - PrintStream's format method can print numbers to two decimal places
   public void printSortedHits(String word) {
      double[] times = new double[numNodes];
      Arrays.fill(times, -1.00);
      double timeAvg = -1.00;
      int count = 0;
   
      for (int i = 0; i < numNodes; i++) {
         for (int j = 1; j < numNodes; j++) {
            if ((adjMatrix[i][j] != null) && (adjMatrix[i][j].getLabel().equals(word))) {
               timeAvg = ((nodeTimes[j] + nodeTimes[i]) / 2);
               times[count] = timeAvg;
               count++;
               }
            }
         }
         
      Arrays.sort(times, 0, count);
      
      for (int i = 0; i < count; i++) {
         if (times[i] != -1.00) {
            System.out.printf("%.2f ", times[i]);
            }
         else {
            System.out.println("Something broked. ");
            }
         }
      
      System.out.print("\n");
      }
   }

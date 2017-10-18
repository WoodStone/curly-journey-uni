import scala.annotation.tailrec
import scala.collection.mutable
import scala.collection.mutable.ListBuffer
import scala.io.Source

object Main extends App {
  if (args.length != 2) {
    println("Usage: scala Main.scala <bfs|dijkstra|astar> <input>")
    System.exit(1)
  }

  val debug: Boolean = false
  val algorithm: String = args(0)
  val inputFile: String = args(1)

  val matrix: Matrix = Matrix(inputFile)
  val (start, end): (Start, End) = matrix.constructTree()


  println("\nINPUT")
  matrix.out()

  algorithm match {
    case "bfs" => GraphFunctions.BFS(start, end)
    case "dijkstra" => GraphFunctions.Dijkstra(start, end)
    case "astar" => GraphFunctions.AStar(start, end, HeuristicFunctions.Manhattan)
  }

  setPathFlag(end)

  println("\nPATH")
  matrix.out()
  println()


  @tailrec def setPathFlag(current: Cell): Unit = {
    current.path = true
    current.cameFrom match {
      case Some(cell) => setPathFlag(cell)
      case None => /**/
    }
  }

  def setOpenFlag(openSet: List[Cell], closedSet: List[Cell]): Unit = {
    if (debug) {
      openSet.foreach(cell => cell.open = Some(true))
      closedSet.foreach(cell => cell.open = Some(false))
    }
  }

  def printOpenClosed(openSet: List[Cell], closedSet: List[Cell]): Unit = {
    println("\nDATA")
    println("Openset: " + openSet.length)
    println("Closedset: " + closedSet.size)
  }

  case class Matrix(resource: String) {
    val (data, xDims, yDims): (Array[Array[Cell]], Int, Int) = generateNodes()

    private def getDims(line: String): (Int, Int) = {
      val dims = line.split(" ")
      (dims(0).toInt, dims(1).toInt)
    }

    private def generateNodes(): (Array[Array[Cell]], Int, Int) = {
//      val lines: Iterator[String] = Source.fromResource(resource).getLines()
      val lines: Iterator[String] = Source.fromFile(resource).getLines()
      val (xDims, yDims) = getDims(lines.next())
      val data = Array.ofDim[Cell](yDims, xDims)

      lines.zipWithIndex.foreach {
        case (line, y) => {
          line.split("").zipWithIndex.foreach {
            case (c, x) => {
              c match {
                case "." => data(y)(x) = Open(x, y, 1)
                case "#" => data(y)(x) = Closed(x, y)
                case "A" => data(y)(x) = Start(x, y)
                case "B" => data(y)(x) = End(x, y)
                case "w" => data(y)(x) = Water(x, y)
                case "m" => data(y)(x) = Mountains(x, y)
                case "f" => data(y)(x) = Forests(x, y)
                case "g" => data(y)(x) = Grasslands(x, y)
                case "r" => data(y)(x) = Roads(x, y)
              }
            }
          }
        }
      }

      (data, xDims, yDims)
    }

    def constructTree(): (Start, End) = {
      //TODO rewrite
      var start: Start = null
      var end: End = null
      for ((line, y) <- data.zipWithIndex) {
        line.zipWithIndex.foreach {
          case (cell, x) => {
            cell match {
              case c: Start => start = c
              case c: End => end = c
              case _ => /**/
            }
            cell.neighbors = findNeighbors(cell)
          }
        }
      }

      (start, end)
    }

    def findNeighbors(cell: Cell): List[Cell] = {
      val list: mutable.ListBuffer[Cell] = new ListBuffer[Cell]
      val (x, y) = (cell.x, cell.y)

      if (x - 1 >= 0) list += data(y)(x - 1)
      if (x + 1 < xDims) list += data(y)(x + 1)
      if (y - 1 >= 0) list += data(y - 1)(x)
      if (y + 1 < yDims) list += data(y + 1)(x)

      list.toList
    }

    def out(): Unit = {
      data.foreach(line => {
        line.foreach {
          case c => {
            if (c.path) {
              print("O")
            } else if (c.open.isDefined) {
              if (c.open.get) print("*") else print("X")
            } else {
              print(c.symbol)
            }
          }
        }
        println()
      })
    }

  }

  object GraphFunctions {

    def AStar(start: Start, end: End, heuristic_func: (Cell, Cell) => Double): Boolean = {
      search(start, end, Some(heuristic_func), OrderingFunctions.LowestFScore)
    }

    def Dijkstra(start: Start, end: End): Boolean = {
      search(start, end, None, OrderingFunctions.LowestGScore)
    }

    def BFS(start: Start, end: End): Boolean = {
      val closedSet: mutable.Set[Cell] = new mutable.HashSet[Cell]()
      val openSet: mutable.Queue[Cell] = mutable.Queue.empty[Cell]

      openSet.enqueue(start)

      while (openSet.nonEmpty) {
        val current: Cell = openSet.dequeue()

        if (current.isInstanceOf[End]) {
          /*Finished*/
          printOpenClosed(openSet.toList, closedSet.toList)
          setOpenFlag(openSet.toList, closedSet.toList)
          return true
        }

        closedSet += current

        current.neighbors.filter(n => {
          if (n.isInstanceOf[Closed]) closedSet += n
          !closedSet.contains(n)
        }).foreach(neighbor => {
          if (!openSet.contains(neighbor)) openSet += neighbor
          neighbor.cameFrom = Some(current)
        })
      }

      return false
    }

    private def search(start: Start, end: End, heuristic_func: Option[(Cell, Cell) => Double], ordering_func: Ordering[Cell]): Boolean = {
      val closedSet: mutable.Set[Cell] = mutable.HashSet.empty[Cell]
      val openSet: java.util.PriorityQueue[Cell] = new java.util.PriorityQueue[Cell](11, ordering_func)

      start.gScore = 0
      heuristic_func.foreach(heuristic => start.fScore = heuristic(start, end))
//      start.fScore = heuristic_func(start, end)
      openSet.add(start)

      while (!openSet.isEmpty) {
        val current: Cell = openSet.remove()

        if (current.isInstanceOf[End]) {
          /*Finished*/
          printOpenClosed(openSet.toArray().toList.asInstanceOf[List[Cell]], closedSet.toList)
          setOpenFlag(openSet.toArray().toList.asInstanceOf[List[Cell]], closedSet.toList)
          return true
        }

        closedSet += current

        current.neighbors.filter(n => {
          if (n.isInstanceOf[Closed]) closedSet += n
          !closedSet.contains(n)
        }).filter(neighbor => {
          var result: Boolean = true
          if (!openSet.contains(neighbor)) {
            openSet.add(neighbor)
          }
          neighbor.tScore = current.gScore + neighbor.cost
          if (neighbor.tScore >= neighbor.gScore && neighbor.gScore != -1) {
            result = false
          }
          result
        }).foreach(neighbor => {
          neighbor.cameFrom = Some(current)
          neighbor.gScore = neighbor.tScore
          heuristic_func.foreach(heuristic => neighbor.fScore = neighbor.gScore + heuristic(neighbor, end))
//          neighbor.fScore = neighbor.gScore + heuristic_func(neighbor, end)
          openSet.remove(neighbor)
          openSet.add(neighbor)

        })
      }

      return false
    }

  }

  object OrderingFunctions {
    def LowestFScore(a: Cell, b: Cell): Int = a.fScore.compare(b.fScore)
    def LowestGScore(a: Cell, b: Cell): Int = a.gScore.compare(b.gScore)
  }

  object HeuristicFunctions {

    def Manhattan(a: Cell, b: Cell): Double = {
      val xDelta: Int = Math.abs(a.x - b.x)
      val yDelta: Int = Math.abs(a.y - b.y)
      xDelta + yDelta
    }

    def Euclidean(a: Cell, b: Cell): Double = {
      val xComp: Double = Math.pow(a.x - b.x, 2)
      val yComp: Double = Math.pow(a.y - b.y, 2)
      Math.sqrt(xComp + yComp)
    }

  }

  class Cell(val x: Int, val y: Int, val symbol: String, val cost: Int = 1) {
    var neighbors: List[Cell] = _
    var cameFrom: Option[Cell] = None
    var fScore: Double = -1
    var gScore: Double = -1
    var tScore: Double = -1
    var path: Boolean = false
    var open: Option[Boolean] = None
  }

  case class Open(override val x: Int, override val y: Int, override val cost: Int) extends Cell(x, y, ".", cost)
  case class Closed(override val x: Int, override val y: Int) extends Cell(x, y, "#")
  case class Start(override val x: Int, override val y: Int) extends Cell(x, y, "A")
  case class End(override val x: Int, override val y: Int) extends Cell(x, y, "B")
  case class Water(override val x: Int, override val y: Int) extends Cell(x, y, "w", 100)
  case class Mountains(override val x: Int, override val y: Int) extends Cell(x, y, "m", 50)
  case class Forests(override val x: Int, override val y: Int) extends Cell(x, y, "f", 10)
  case class Grasslands(override val x: Int, override val y: Int) extends Cell(x, y, "g", 5)
  case class Roads(override val x: Int, override val y: Int) extends Cell(x, y, "r", 1)

}

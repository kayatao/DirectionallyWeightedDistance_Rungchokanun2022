/* Exemplary file for Chapter 6 - Exploring Graphs. */

using System.Collections.Generic;

namespace Graphs
{
    public class Node<T>
    {
        public int Index { get; set; }
        public T Data { get; set; }
        public List<Node<T>> Neighbors { get; set; } = new List<Node<T>>();
        public List<int> Weights { get; set; } = new List<int>();

        public int level;

        public bool visited;

        public string label;
        public List<double> r { get; set; } = new List<double>();
        public List<double> theta { get; set; } = new List<double>();

        public override string ToString()
        {
            return $"Node with index {Index}: {Data}, label: {label}, neighbors: {Neighbors.Count}";
        }
    }
}

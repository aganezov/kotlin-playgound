package ensembl

import java.io.BufferedReader
import java.io.File
import java.io.FileReader

val ID_KEY = "id"
val START_KEY = "start"
val END_KEY = "end"
val STRAND_KEY = "strand"
val CHROMOSOME_NAME_KEY = "chr"

val headerMapper = mapOf(
        // Gene id headers list
        "Ensembl Gene ID" to ID_KEY,

        // Gene start headers list
        "Gene Start (bp)" to START_KEY,

        // Gene end headers list
        "Gene End (bp)" to END_KEY,

        // Gene strand headers list
        "Strand" to STRAND_KEY,

        // Chromosome name headers list
        "Chromosome Name" to CHROMOSOME_NAME_KEY
)

fun transformHeaders(headers: Iterable<String>) = headers.map { headerMapper.getOrElse(it, { "-1" }) }

val FILE_PATH = "/Users/aganezov/Desktop"
val FILE_NAME = "chr14_genes_ensembl.txt"
val SEPARATOR = "\t"

fun main(args: Array<String>) {
    val genes = readGeneOrder(filePath = FILE_PATH, fileName = FILE_NAME, separator = SEPARATOR).sortedBy { it.start }
    val overlappingGenes = genes.dropLast(1).zip(genes.drop(1)).filter { it.first.end > it.second.start }
    println(overlappingGenes.size)

}

fun readGeneOrder(filePath: String, fileName: String, separator: String = "\t"): List<Gene> {
    val reader = BufferedReader(FileReader(filePath + File.separator + fileName))
    val headers = transformHeaders(reader.readLine().trim().split(separator))
    return reader.lineSequence().map { fromCSVString(CSVString = it, headers = headers, separator = separator) }.toList()
}

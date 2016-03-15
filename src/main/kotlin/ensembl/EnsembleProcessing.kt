package ensembl

import java.io.BufferedReader
import java.io.File
import java.io.FileReader

val headerMapper = mapOf<String, String>(
        // Gene id headers list
        "Ensembl Gene ID" to "id",

        // Gene start headers list
        "Gene Start (bp)" to "start",

        // Gene end headers list
        "Gene End (bp)" to "end",

        // Gene strand headers list
        "Strand" to "strand",

        // Chromosome name headers list
        "Chromosome Name" to "chr"
)

fun transformHeaders(headers: Iterable<String>) = headers.map { headerMapper.getOrElse(it, { "-1" }) }

val FILE_PATH = "/Users/aganezov/Desktop"
val FILE_NAME = "chr14_genes_ensembl.txt"
val SEPARATOR = "\t"

fun main(args: Array<String>) {
    val genes = readGeneOrder(filePath = FILE_PATH, fileName = FILE_NAME, separator = SEPARATOR)
            .sortedBy { it.start }
    val overlappingGenes = genes.dropLast(1).zip(genes.drop(1)).filter { it.first.end > it.second.start }

}

fun readGeneOrder(filePath: String, fileName: String, separator: String = "\t"): List<Gene> {
    val reader = BufferedReader(FileReader(filePath + File.separator + fileName))
    val headers = transformHeaders(reader.readLine().trim().split(separator))
    return reader.lineSequence().map { fromCSVString(CSVString = it, headers = headers, separator = separator) }.toList()
}

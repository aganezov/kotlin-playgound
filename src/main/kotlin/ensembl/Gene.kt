package ensembl

fun fromCSVString(CSVString: String, headers: Iterable<String>, separator: String = "\t"): Gene {
    val data = headers.zip(CSVString.split(separator)).map { it.first to it.second }.toMap()
    return Gene(
            id = data["id"] ?: "-1",
            start = data["start"]?.toLong() ?: -1L,
            end = data["end"]?.toLong() ?: -1L,
            strand = data["strand"]?.toInt() ?: 0,
            chromosomeName = data["chr"] ?: "-1"
    )
}

data class Gene(
        val id: String,
        val start: Long,
        val end: Long,
        val strand: Int,
        val chromosomeName: String
)


